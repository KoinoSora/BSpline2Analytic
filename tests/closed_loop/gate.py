"""Deterministic evidence rules; no SDK import, build, subprocess or network access.

The caller must supply policy, registry and live context from its trusted control
checkout. Hashes detect corruption, not an untrusted producer forging results.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
import re
import statistics


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def positive(value):
    return type(value) in (float, int) and math.isfinite(value) and value > 0


def ratio_interval(before, after, policy):
    """Paired bootstrap of the ratio of medians, with reproducible resampling.

    Pair i must be the same alternating A/B measurement block. Constant fixtures
    give an exact interval. This does not model scheduling or thermal drift.
    """
    point = statistics.median(after) / statistics.median(before)
    rng = random.Random(policy["bootstrap_seed"])
    n = len(before)
    ratios = []
    for _ in range(policy["bootstrap_samples"]):
        sample = [rng.randrange(n) for _ in range(n)]
        ratios.append(statistics.median(after[i] for i in sample) /
                      statistics.median(before[i] for i in sample))
    ratios.sort()
    tail = (1 - policy["confidence"]) / 2
    return point, ratios[int(tail * (len(ratios) - 1))], ratios[math.ceil((1-tail) * (len(ratios)-1))]


def evaluate(policy, registry, evidence, context):
    """Never return production PASS for simulated evidence or an inactive policy."""
    issues, warnings, comparisons = [], [], []

    def fail(code, detail):
        issues.append({"code": code, "detail": detail})

    simulation = evidence.get("mode") == "simulation"
    if evidence.get("mode") not in ("simulation", "production"):
        fail("MODE", "Unknown evidence mode")
    if not simulation:
        if not policy.get("real_tests_enabled"):
            fail("NOT_ACTIVATED", "Real tests and the production gate are disabled")
        if registry.get("state") != "APPROVED" or not registry.get("approval"):
            fail("BASELINE_NOT_APPROVED", "No administrator-approved baseline")
        if context.get("trusted_collector") is not True:
            fail("UNTRUSTED_COLLECTOR", "Producer identity has not been verified by the controller")
        if evidence.get("synthetic") is not False:
            fail("SYNTHETIC_PRODUCTION", "Synthetic or unclassified data cannot authorize a merge")
    elif evidence.get("synthetic") is not True:
        fail("SIMULATION_MARKER", "Simulation must be explicitly marked synthetic")
    for key in ("head_sha", "base_sha", "merge_tree_sha"):
        if not re.fullmatch(r"[0-9a-f]{40}", str(context.get(key, ""))):
            fail("CONTEXT", f"Invalid {key}")
        if evidence.get(key) != context.get(key):
            fail("STALE_EVIDENCE", f"{key} no longer matches the live PR")
    if context.get("repository") != policy["repository"]:
        fail("REPOSITORY", "Repository mismatch")
    if context.get("base_branch") != policy["target_branch"]:
        fail("BASE_BRANCH", "PR does not target the protected branch")
    for key, expected in (("policy_sha256", digest(policy)),
                          ("registry_sha256", digest(registry)),
                          ("environment_sha256", registry.get("environment_sha256"))):
        if not expected or evidence.get(key) != expected:
            fail("IDENTITY", f"Mismatched {key}")
    if evidence.get("anchor_sha") != registry.get("anchor_sha"):
        fail("ANCHOR", "Fixed baseline commit changed")
    if evidence.get("timing_scope") != "complete-api-cold-context":
        fail("TIMING_SCOPE", "Primary measurements must cover the full API with fresh context")
    if evidence.get("environment_valid") is not True:
        fail("ENVIRONMENT", "Environment drift or load invalidated this batch")
    docs = evidence.get("documents", {})
    for kind in ("algorithm", "impact"):
        if not isinstance(docs.get(kind), str) or not docs[kind].strip():
            fail("DOCUMENTS", f"Missing {kind} document")

    cases = registry.get("cases", [])
    if not cases:
        fail("EMPTY_REGISTRY", "An empty suite is not a passing suite")
    keys = [x["id"] for x in cases]
    if len(set(keys)) != len(keys):
        fail("DUPLICATE_CASE", "Registry contains duplicate case keys")
    measurements = evidence.get("cases", [])
    measured_keys = [x["id"] for x in measurements]
    if len(set(measured_keys)) != len(measured_keys) or set(measured_keys) != set(keys):
        fail("COVERAGE", "Missing, additional or duplicate measurements")
    measured = {x["id"]: x for x in measurements}
    groups = {}
    for spec in cases:
        key = spec["id"]
        row = measured.get(key)
        if row is None:
            continue
        contract = spec.get("contract")
        if contract not in policy["contracts"]:
            fail("CONTRACT", key)
            continue
        if spec.get("entrypoint") in policy["retired_entrypoints"] or any(
                value in spec.get("dataset", "") for value in policy["retired_datasets"]):
            fail("RETIRED_TEST", key)
        for name in ("input_sha256", "contract", "entrypoint", "options_sha256"):
            if row.get(name) != spec.get(name) or not row.get(name):
                fail("CASE_IDENTITY", f"{key}: {name}")
        runs = row.get("runs", [])
        count = policy["required_repetitions"]
        if len(runs) != count or [x.get("iteration") for x in runs] != list(range(count)):
            fail("REPETITIONS", key)
            continue
        case_times = {name: [] for name in ("candidate_ms", "anchor_ms", "base_ms")}
        valid_timing = True
        for run in runs:
            status = run.get("status")
            expected = spec["expectation"]
            if expected == "diagnostic":
                allowed_diagnostics = {"ACCEPT", "REJECT", "UNRESOLVED", "TIMEOUT", "INVALID", "UNSUPPORTED", "BUDGET_EXCEEDED"}
                if spec.get("diagnostic_status") not in allowed_diagnostics or status != spec.get("diagnostic_status"):
                    fail("DIAGNOSTIC", key)
            elif status in ("TIMEOUT", "CRASH", "NOT_RUN", "INVALID"):
                fail("EXECUTION", f"{key}: {status}")
            elif expected == "accept" and status != "ACCEPT":
                fail("CORRECTNESS", f"{key}: positive not accepted")
            elif expected == "reject" and status != "REJECT":
                fail("FALSE_ACCEPT", f"{key}: negative not rejected")
            elif expected == "known_unresolved" and status not in ("UNRESOLVED", "ACCEPT"):
                fail("DEBT_REGRESSION", key)
            elif expected not in ("accept", "reject", "known_unresolved", "diagnostic"):
                fail("EXPECTATION", key)
            if expected == "reject" and status == "REJECT" and run.get("rejection_evidence") != "PASS":
                fail("REJECTION_EVIDENCE", key)
            if expected == "known_unresolved" and not spec.get("debt_issue"):
                fail("DEBT_RECORD", key)
            if status == "ACCEPT":
                if run.get("independent_audit") != "PASS":
                    fail("AUDIT", key)
                checks = run.get("checks", {})
                if any(checks.get(c) != "PASS" for c in policy["contracts"][contract]):
                    fail("GEOMETRY", f"{key}: missing/failed contract evidence")
                # Scalar bounds belong to this contract; sampling is never labelled continuous.
                if contract == "mother-distance-joint-v1":
                    bounds = run.get("bounds", {})
                    for field, limit in (("surface_distance", 1e-5), ("boundary_position", 1e-5),
                                         ("joint_tangent_deg", 0.1)):
                        value = bounds.get(field)
                        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value < limit:
                            fail("GEOMETRY_BOUND", f"{key}: {field}")
            if spec.get("new") and expected == "known_unresolved":
                fail("NEW_DEBT", key)
            if expected == "diagnostic":
                continue
            for name in case_times:
                value = run.get(name)
                if spec.get("new") and name != "candidate_ms":
                    continue
                if not positive(value):
                    fail("TIMING", f"{key}: missing/nonpositive/nonfinite {name}")
                    valid_timing = False
                else:
                    case_times[name].append(value)
        if spec["expectation"] == "diagnostic" or not valid_timing:
            continue
        after = case_times["candidate_ms"]
        if spec.get("new"):
            if max(after) >= policy["new_case_limit_ms"]:
                fail("NEW_CASE_10MS", key)
            continue
        if policy["preserve_existing_10ms"] and max(case_times["anchor_ms"]) < 10 and max(after) >= 10:
            fail("EXISTING_CASE_10MS", key)
        for baseline in ("anchor_ms", "base_ms"):
            before = case_times[baseline]
            point, lower, upper = ratio_interval(before, after, policy)
            comparisons.append(dict(case=key, baseline=baseline, ratio=point, lower=lower, upper=upper))
            if lower > policy["maximum_time_ratio"]:
                fail("PERFORMANCE", f"{key} vs {baseline}: {point:.6f}")
            elif upper > policy["maximum_time_ratio"]:
                warnings.append({"code": "UNCERTAIN_PERFORMANCE", "detail": f"{key} vs {baseline}"})
            for group in ("all", spec.get("group", "ungrouped")):
                pair = groups.setdefault((group, baseline), [0.0, 0.0])
                pair[0] += statistics.median(before)
                pair[1] += statistics.median(after)
    for (group, baseline), (before, after) in groups.items():
        comparisons.append(dict(group=group, baseline=baseline, ratio=after/before))
        if after/before > policy["maximum_time_ratio"]:
            fail("GROUP_PERFORMANCE", f"{group} vs {baseline}")
    state = "FAIL" if issues else "INCONCLUSIVE" if warnings else "PASS"
    if simulation:
        state = "SIMULATION_" + state
    return dict(schema_version=1, state=state, mode=evidence.get("mode"),
                synthetic=simulation, production_eligible=state == "PASS",
                issues=issues, warnings=warnings, comparisons=comparisons,
                head_sha=context.get("head_sha"), base_sha=context.get("base_sha"),
                merge_tree_sha=context.get("merge_tree_sha"), evidence_sha256=digest(evidence),
                policy_sha256=digest(policy), registry_sha256=digest(registry))


def merge_decision(result, context, publication, reviews, policy, simulation=False):
    """Controller-only final decision; no API that performs a merge is provided."""
    prefix = "SIMULATION_" if simulation else ""
    if result["state"] != prefix + "PASS" or result.get("synthetic") != simulation:
        return prefix + "BLOCKED_VALIDATION"
    if not simulation and (not policy.get("real_tests_enabled") or not result.get("production_eligible")):
        return "NOT_ACTIVATED"
    for key in ("head_sha", "base_sha", "merge_tree_sha"):
        if result.get(key) != context.get(key):
            return prefix + "STALE"
    if publication.get("verified") is not True or publication.get("evidence_sha256") != result["evidence_sha256"]:
        return prefix + "BLOCKED_PUBLICATION"
    latest = {}
    for review in sorted(reviews, key=lambda r: r["sequence"]):
        if review["state"] in ("APPROVED", "CHANGES_REQUESTED", "DISMISSED"):
            latest[review["login"]] = review
    if any(r["state"] == "CHANGES_REQUESTED" for r in latest.values()):
        return prefix + "CHANGES_REQUESTED"
    for login, review in latest.items():
        if (login in policy["reviewers"] and login != context.get("author") and
                review["state"] == "APPROVED" and review["head_sha"] == context["head_sha"]):
            return prefix + "READY_FOR_ADMIN_MERGE"
    return prefix + "AWAITING_ADMIN_REVIEW"
