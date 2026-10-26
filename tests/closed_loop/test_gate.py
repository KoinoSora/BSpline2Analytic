"""Protocol fault injection only. No geometric inputs or SDK are loaded."""
import copy
import json
from pathlib import Path
import tempfile
import unittest

from gate import digest, evaluate, merge_decision
from run import ROOT, fixture, read, smoke, verify


class GateTests(unittest.TestCase):
    def setUp(self):
        self.policy = read(ROOT / "policy.json")
        self.context = dict(repository=self.policy["repository"], base_branch="main",
                            head_sha="1" * 40, base_sha="2" * 40, merge_tree_sha="3" * 40,
                            author="contributor", trusted_collector=False)
        self.registry, self.evidence = fixture(self.policy, self.context)

    def evaluate(self):
        return evaluate(self.policy, self.registry, self.evidence, self.context)

    def assert_code(self, code):
        result = self.evaluate()
        self.assertIn(code, [x["code"] for x in result["issues"]])
        self.assertFalse(result["production_eligible"])

    def update_registry(self):
        self.evidence["registry_sha256"] = digest(self.registry)

    def test_simulation_never_production_pass(self):
        self.assertEqual(self.evaluate()["state"], "SIMULATION_PASS")
        self.assertFalse(self.evaluate()["production_eligible"])

    def test_real_mode_disabled(self):
        self.evidence.update(mode="production", synthetic=False)
        for code in ("NOT_ACTIVATED", "BASELINE_NOT_APPROVED", "UNTRUSTED_COLLECTOR"):
            self.assert_code(code)

    def test_synthetic_cannot_claim_production(self):
        self.evidence["mode"] = "production"
        self.assert_code("SYNTHETIC_PRODUCTION")

    def test_slow_existing_case(self):
        for run in self.evidence["cases"][0]["runs"]:
            run["candidate_ms"] = 4.24
        self.assert_code("PERFORMANCE")

    def test_exact_five_percent_boundary(self):
        for run in self.evidence["cases"][0]["runs"]:
            run["candidate_ms"] = 4.2
        self.assertEqual(self.evaluate()["state"], "SIMULATION_PASS")

    def test_anchor_prevents_cumulative_drift(self):
        for run in self.evidence["cases"][0]["runs"]:
            run.update(candidate_ms=4.3, base_ms=4.2)
        self.assert_code("PERFORMANCE")

    def test_fast_average_does_not_hide_single_regression(self):
        other = copy.deepcopy(self.registry["cases"][0])
        other["id"] = "synthetic:fast-other"
        row = copy.deepcopy(self.evidence["cases"][0])
        row["id"] = other["id"]
        for run in row["runs"]:
            run["candidate_ms"] = 1
        self.registry["cases"].append(other)
        self.evidence["cases"].append(row)
        self.update_registry()
        for run in self.evidence["cases"][0]["runs"]:
            run["candidate_ms"] = 4.4
        self.assert_code("PERFORMANCE")

    def test_new_case_exactly_ten_fails(self):
        self.evidence["cases"][1]["runs"][0]["candidate_ms"] = 10
        self.assert_code("NEW_CASE_10MS")

    def test_original_ten_ms_boundary(self):
        for run in self.evidence["cases"][0]["runs"]:
            run.update(anchor_ms=9.8, base_ms=9.8, candidate_ms=10.1)
        self.assert_code("EXISTING_CASE_10MS")

    def test_missing_case(self):
        self.evidence["cases"].pop()
        self.assert_code("COVERAGE")

    def test_duplicate_case(self):
        self.evidence["cases"].append(self.evidence["cases"][0])
        self.assert_code("COVERAGE")

    def test_empty_suite(self):
        self.registry["cases"] = []
        self.evidence["cases"] = []
        self.update_registry()
        self.assert_code("EMPTY_REGISTRY")

    def test_incomplete_repetitions(self):
        self.evidence["cases"][0]["runs"].pop()
        self.assert_code("REPETITIONS")

    def test_nonfinite_or_zero_times(self):
        for number in (float("nan"), float("inf"), 0, -1, True):
            with self.subTest(value=number):
                self.evidence["cases"][0]["runs"][0]["candidate_ms"] = number
                # Canonical hashing also rejects NaN/Inf instead of sealing it.
                if number != number or number == float("inf"):
                    with self.assertRaises(ValueError):
                        self.evaluate()
                else:
                    self.assert_code("TIMING")

    def test_false_accept(self):
        self.registry["cases"][0]["expectation"] = "reject"
        self.update_registry()
        self.assert_code("FALSE_ACCEPT")

    def test_timeout_is_not_correct_rejection(self):
        self.registry["cases"][0]["expectation"] = "reject"
        self.update_registry()
        self.evidence["cases"][0]["runs"][0]["status"] = "TIMEOUT"
        self.assert_code("EXECUTION")

    def test_negative_needs_rejection_evidence(self):
        self.registry["cases"][0]["expectation"] = "reject"
        self.update_registry()
        for run in self.evidence["cases"][0]["runs"]:
            run["status"] = "REJECT"
        self.assert_code("REJECTION_EVIDENCE")
        for run in self.evidence["cases"][0]["runs"]:
            run["rejection_evidence"] = "PASS"
        self.assertEqual(self.evaluate()["state"], "SIMULATION_PASS")

    def test_intentional_timeout_diagnostic_is_separate(self):
        self.registry["cases"][0].update(expectation="diagnostic", diagnostic_status="TIMEOUT")
        self.update_registry()
        for run in self.evidence["cases"][0]["runs"]:
            run.update(status="TIMEOUT", candidate_ms=None, base_ms=None, anchor_ms=None)
        result = self.evaluate()
        self.assertEqual(result["state"], "SIMULATION_PASS")
        self.assertEqual(result["comparisons"], [])

    def test_missing_diagnostic_cannot_be_labelled_expected(self):
        for status in (None, "NOT_RUN", "CRASH"):
            with self.subTest(status=status):
                self.registry["cases"][0].update(expectation="diagnostic", diagnostic_status=status)
                self.update_registry()
                for run in self.evidence["cases"][0]["runs"]:
                    run["status"] = status
                self.assert_code("DIAGNOSTIC")

    def test_geometric_bound_cannot_be_relabelled_pass(self):
        self.evidence["cases"][0]["runs"][0]["bounds"]["surface_distance"] = 1e-5
        self.assert_code("GEOMETRY_BOUND")

    def test_missing_independent_audit(self):
        self.evidence["cases"][0]["runs"][0]["independent_audit"] = "NOT_RUN"
        self.assert_code("AUDIT")

    def test_missing_contract_check(self):
        del self.evidence["cases"][0]["runs"][0]["checks"]["uv_joint_tangent"]
        self.assert_code("GEOMETRY")

    def test_mismatched_contract(self):
        self.evidence["cases"][0]["contract"] = "general-rational-v2"
        self.assert_code("CASE_IDENTITY")

    def test_unresolved_cannot_be_new_success(self):
        self.registry["cases"][1].update(expectation="known_unresolved", debt_issue="synthetic-issue")
        self.update_registry()
        self.assert_code("NEW_DEBT")

    def test_known_unresolved_recovery_still_subject_to_performance(self):
        self.registry["cases"][0].update(expectation="known_unresolved", debt_issue="synthetic-issue")
        self.update_registry()
        for run in self.evidence["cases"][0]["runs"]:
            run["candidate_ms"] = 8
        self.assert_code("PERFORMANCE")

    def test_retired_entrypoint(self):
        self.registry["cases"][0]["entrypoint"] = "ValidateGeneralHybrid"
        self.update_registry()
        self.assert_code("RETIRED_TEST")

    def test_retired_dataset(self):
        self.registry["cases"][0]["dataset"] = "tests/data/general-extreme320/x"
        self.update_registry()
        self.assert_code("RETIRED_TEST")

    def test_stale_source_or_base(self):
        for key in ("head_sha", "base_sha", "merge_tree_sha"):
            with self.subTest(key=key):
                previous = self.context[key]
                self.context[key] = "f" * 40
                self.assert_code("STALE_EVIDENCE")
                self.context[key] = previous

    def test_modified_policy(self):
        self.policy["maximum_time_ratio"] = 1.5
        self.assert_code("IDENTITY")

    def test_wrong_environment(self):
        self.evidence["environment_sha256"] = "bad"
        self.assert_code("IDENTITY")

    def test_missing_document(self):
        self.evidence["documents"].pop("impact")
        self.assert_code("DOCUMENTS")

    def test_noise_is_inconclusive(self):
        for i, run in enumerate(self.evidence["cases"][0]["runs"]):
            run["candidate_ms"] = 4.0 if i < 11 else 4.4
        self.assertEqual(self.evaluate()["state"], "SIMULATION_INCONCLUSIVE")

    def test_review_lifecycle(self):
        result = self.evaluate()
        publication = dict(verified=True, evidence_sha256=result["evidence_sha256"])
        approval = dict(sequence=1, login=self.policy["reviewers"][0], state="APPROVED",
                        head_sha=self.context["head_sha"])
        decision = lambda reviews: merge_decision(result, self.context, publication, reviews, self.policy, True)
        self.assertEqual(decision([]), "SIMULATION_AWAITING_ADMIN_REVIEW")
        self.assertEqual(decision([approval]), "SIMULATION_READY_FOR_ADMIN_MERGE")
        self.assertEqual(decision([dict(approval, login="not-an-admin")]), "SIMULATION_AWAITING_ADMIN_REVIEW")
        self.assertEqual(decision([dict(approval, head_sha="a" * 40)]), "SIMULATION_AWAITING_ADMIN_REVIEW")
        self.assertEqual(decision([approval, dict(approval, sequence=2, state="CHANGES_REQUESTED")]),
                         "SIMULATION_CHANGES_REQUESTED")
        self.assertEqual(decision([approval, dict(approval, sequence=2, state="DISMISSED")]),
                         "SIMULATION_AWAITING_ADMIN_REVIEW")
        self.context["author"] = approval["login"]
        self.assertEqual(decision([approval]), "SIMULATION_AWAITING_ADMIN_REVIEW")
        self.assertEqual(merge_decision(result, self.context, publication, [approval], self.policy), "BLOCKED_VALIDATION")

    def test_full_offline_lifecycle_and_tamper(self):
        with tempfile.TemporaryDirectory() as directory:
            bundle = Path(directory) / "bundle"
            smoke(bundle, self.context, self.policy)
            self.assertEqual(len(verify(bundle)["files"]), 9)
            (bundle / "report.md").write_text("modified", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify(bundle)


if __name__ == "__main__":
    unittest.main()
