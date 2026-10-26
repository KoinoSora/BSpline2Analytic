"""Offline protocol rehearsal and evidence interchange. Never calls SGK."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import sys

from gate import digest, evaluate, merge_decision

ROOT = Path(__file__).resolve().parent


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def fixture(policy, context):
    """These small invented numbers are protocol inputs, NEVER benchmark data."""
    spec = dict(id="synthetic:existing-cone", input_sha256=digest("invented-input"),
                options_sha256=digest("invented-options"), contract="mother-distance-joint-v1",
                entrypoint="ValidateMotherDistanceConversion", dataset="synthetic-only",
                expectation="accept", new=False, group="cone")
    new = dict(spec, id="synthetic:new-cone", new=True)
    registry = dict(schema_version=1, state="SIMULATION_ONLY", approval=None,
                    anchor_sha="a" * 40, environment_sha256=digest("invented-environment"), cases=[spec, new])
    evidence = dict(schema_version=1, mode="simulation", synthetic=True,
                    anchor_sha=registry["anchor_sha"], policy_sha256=digest(policy),
                    registry_sha256=digest(registry), environment_sha256=registry["environment_sha256"],
                    environment_valid=True, timing_scope="complete-api-cold-context",
                    documents={"algorithm": "docs/algorithm/CLOSED_LOOP_VALIDATION.md",
                               "impact": "docs/requirements/CLOSED_LOOP_VALIDATION.md"},
                    cases=[], **{k: context[k] for k in ("head_sha", "base_sha", "merge_tree_sha")})
    for item in registry["cases"]:
        row = {k: item[k] for k in ("id", "input_sha256", "options_sha256", "contract", "entrypoint")}
        row["runs"] = [dict(iteration=i, status="ACCEPT", independent_audit="PASS",
                            checks={k: "PASS" for k in policy["contracts"][item["contract"]]},
                            bounds=dict(surface_distance=2e-6, boundary_position=3e-6, joint_tangent_deg=0.04),
                            candidate_ms=6.0 if item["new"] else 4.04, anchor_ms=4.0, base_ms=4.0)
                       for i in range(policy["required_repetitions"])]
        evidence["cases"].append(row)
    return registry, evidence


def seal(folder):
    folder = Path(folder)
    files = {p.relative_to(folder).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted(folder.rglob("*")) if p.is_file() and p.name != "manifest.json"}
    write(folder / "manifest.json", {"schema_version": 1, "synthetic": True, "files": files})


def verify(folder):
    folder = Path(folder).resolve()
    manifest = read(folder / "manifest.json")
    actual = {p.relative_to(folder).as_posix() for p in folder.rglob("*")
              if p.is_file() and p != folder / "manifest.json"}
    if actual != set(manifest["files"]):
        raise ValueError("Bundle file list changed")
    for name, expected in manifest["files"].items():
        original = folder / name
        path = original.resolve()
        if not path.is_relative_to(folder) or original.is_symlink():
            raise ValueError("Unsafe bundle path")
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError("Bundle hash mismatch: " + name)
    return manifest


def smoke(folder, context, policy=None):
    policy = policy or read(ROOT / "policy.json")
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    registry, evidence = fixture(policy, context)
    result = evaluate(policy, registry, evidence, context)
    if result["state"] != "SIMULATION_PASS" or result["production_eligible"]:
        raise ValueError("Valid simulation failed or escaped into production")
    publication = dict(verified=True, evidence_sha256=result["evidence_sha256"])
    # The author identity in this synthetic lifecycle is deliberately invented.
    review_context = dict(context, author="synthetic-contributor")
    admin = policy["reviewers"][0]
    review = dict(sequence=1, login=admin, state="APPROVED", head_sha=context["head_sha"])
    lifecycle = [dict(step="evaluated", state=result["state"]),
                 dict(step="not_uploaded", state=merge_decision(result, review_context, {}, [], policy, True)),
                 dict(step="published_simulation", state=merge_decision(result, review_context, publication, [], policy, True)),
                 dict(step="simulated_admin_approval", state=merge_decision(result, review_context, publication, [review], policy, True))]
    expected = ["SIMULATION_PASS", "SIMULATION_BLOCKED_PUBLICATION",
                "SIMULATION_AWAITING_ADMIN_REVIEW", "SIMULATION_READY_FOR_ADMIN_MERGE"]
    if [x["state"] for x in lifecycle] != expected:
        raise ValueError("Lifecycle transition failed")
    stale = dict(review_context, base_sha="b" * 40)
    if stale["base_sha"] == context["base_sha"]:
        stale["base_sha"] = "c" * 40
    lifecycle.append(dict(step="base_changed", state=merge_decision(result, stale, publication, [review], policy, True)))
    if lifecycle[-1]["state"] != "SIMULATION_STALE":
        raise ValueError("Stale validation was accepted")
    production = copy.deepcopy(evidence)
    production.update(mode="production", synthetic=False)
    blocked = evaluate(policy, registry, production, context)
    if not any(x["code"] == "NOT_ACTIVATED" for x in blocked["issues"]):
        raise ValueError("Production activation guard failed")
    for name, value in (("context.json", context), ("policy.json", policy),
                        ("registry.json", registry), ("evidence.json", evidence),
                        ("result.json", result), ("lifecycle.json", lifecycle),
                        ("production-refusal.json", blocked)):
        write(folder / name, value)
    with (folder / "measurements.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["synthetic", "case", "iteration", "candidate_ms", "anchor_ms", "base_ms"])
        for row in evidence["cases"]:
            for run in row["runs"]:
                writer.writerow([True, row["id"], run["iteration"], run["candidate_ms"], run["anchor_ms"], run["base_ms"]])
    (folder / "report.md").write_text(
        "# 闭环协议演练（全部为合成数据）\n\n"
        "**没有运行 SGK、几何数据集、性能基准或模型 AGENT；不能用于批准实际代码合并。**\n\n"
        f"PR head: `{context['head_sha']}`\n\nBase: `{context['base_sha']}`\n\n"
        f"Merge tree: `{context['merge_tree_sha']}`\n\n"
        "演练覆盖：规则判断、报告和哈希封装、待发布阻断、等待审批、模拟审批、主分支变化失效。\n"
        "管理员审批是本地状态机模拟；真实管理员尚未批准。远端发布以单独回执为准。\n\n"
        "正式策略仍为关闭，未批准基线；production-refusal.json 记录拒绝原因。\n",
        encoding="utf-8")
    seal(folder)
    verify(folder)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    command = sub.add_parser("smoke", help="Synthetic protocol only; no benchmark execution")
    command.add_argument("--output", type=Path, required=True)
    command.add_argument("--head", required=True)
    command.add_argument("--base", required=True)
    command.add_argument("--merge-tree", required=True)
    command.add_argument("--repository", default="KoinoSora/BSpline2Analytic")
    command = sub.add_parser("verify")
    command.add_argument("--bundle", type=Path, required=True)
    command = sub.add_parser("evaluate")
    for name in ("policy", "registry", "evidence", "context", "output"):
        command.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "smoke":
            context = dict(repository=args.repository, base_branch="main", head_sha=args.head,
                           base_sha=args.base, merge_tree_sha=args.merge_tree, trusted_collector=False,
                           author="synthetic-contributor")
            print(smoke(args.output, context)["state"])
        elif args.command == "verify":
            print("BUNDLE_VERIFIED", len(verify(args.bundle)["files"]))
        else:
            result = evaluate(read(args.policy), read(args.registry), read(args.evidence), read(args.context))
            write(args.output, result)
            print(result["state"])
            return 0 if result["state"] in ("PASS", "SIMULATION_PASS") else 2
        return 0
    except (ValueError, TypeError, KeyError, OSError) as error:
        print("PROTOCOL_ERROR:", error, file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
