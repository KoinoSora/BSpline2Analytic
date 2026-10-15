"""Audit and summarize paired/third-validator sphere-torus measurements (stdlib only)."""
import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median

STAGES = ["preparation_ms", "position_ms", "subdivision_ms", "boundary_ms", "reporting_ms"]
METRICS = STAGES + ["classification_ms", "validation_ms", "algorithm_ms", "source_extraction_ms",
    "source_boundary_ms", "target_preparation_ms", "mapping_ms", "sample_ms", "residual_build_ms",
    "residual_bound_ms", "boundary_position_ms", "boundary_angle_ms"]


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = list(csv.DictReader(args.csv.open(encoding="utf-8-sig")))
    keys = {(int(r["id"]), int(r["run"]), r["backend"]) for r in rows}
    expected = {(i, run, backend) for i in range(301, 500, 2) for run in range(3)
        for backend in ("paired", "mother")}
    assert len(rows) == 600 and keys == expected, "Missing/duplicate test measurements"
    groups, cases = defaultdict(list), defaultdict(list)
    for r in rows:
        values = {k: float(r[k]) for k in METRICS}
        assert all(math.isfinite(v) and v >= 0 for v in values.values())
        assert abs(sum(values[k] for k in STAGES) - values["validation_ms"]) < 1e-7
        assert abs(values["algorithm_ms"] - values["classification_ms"] - values["validation_ms"]) < 1e-7
        if r["status"] != "Passed":
            assert int(r["collapsed"]) and r["status"] in ("NoCandidate", "ClassificationFailed")
            assert values["validation_ms"] == 0
            continue
        assert not int(r["collapsed"]) and int(r["source_cache_hits"]) == 0
        assert float(r["mother_or_paired_bound"]) < 1e-5
        assert float(r["tangent_bound"]) < .1 and float(r["normal_bound"]) < .1
        if r["backend"] == "mother":
            assert float(r["boundary_position_bound"]) < 1e-5
            assert sum(values[k] for k in ("sample_ms", "residual_build_ms", "residual_bound_ms")) <= values["position_ms"] + 1e-5
            assert values["boundary_position_ms"] + values["boundary_angle_ms"] <= values["boundary_ms"] + 1e-5
        groups[(r["type"], r["backend"])].append(values)
        cases[(int(r["id"]), r["type"], r["backend"])].append(values)
    per_case = []
    for (i, surface, backend), runs in sorted(cases.items()):
        assert len(runs) == 3
        per_case.append(dict(id=i, type=surface, backend=backend, runs=3,
            **{k: mean(r[k] for r in runs) for k in METRICS},
            algorithm_median_ms=median(r["algorithm_ms"] for r in runs),
            algorithm_max_run_ms=max(r["algorithm_ms"] for r in runs)))
    modules = []
    for (surface, backend), runs in sorted(groups.items()):
        case_rows = [r for r in per_case if r["type"] == surface and r["backend"] == backend]
        averages = {k: mean(r[k] for r in runs) for k in METRICS}
        alg = sorted(r["algorithm_ms"] for r in case_rows)
        modules.append(dict(type=surface, backend=backend, cases=len(case_rows), measurements=len(runs),
            **averages, **{k.replace("_ms", "_percent"): 100 * averages[k] / averages["validation_ms"] for k in STAGES},
            algorithm_p95_case_mean_ms=alg[math.ceil(.95 * len(alg)) - 1],
            algorithm_max_case_mean_ms=max(alg),
            cases_mean_under_10ms=sum(x < 10 for x in alg),
            cases_all_runs_under_10ms=sum(r["algorithm_max_run_ms"] < 10 for r in case_rows)))
    summary = dict(schema=1, input_rows=len(rows), runs_per_case_backend=3,
        contract="mother distance with floating-point guard plus interval paired boundaries",
        algorithm_ms_scope="classification compute plus first candidate validation; no optimizer run",
        status_measurements=dict(Counter(f'{r["type"]}/{r["backend"]}/{r["status"]}' for r in rows)),
        unique_no_candidate_ids=sorted({int(r["id"]) for r in rows if r["status"] != "Passed"}),
        modules=modules)
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / "modules-by-type.csv", modules)
    write_csv(args.output / "per-case.csv", per_case)
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
