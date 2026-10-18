# 文件职责：母面/配对验证记录的分组统计。
# 识别历史与当前角度合同，不能将联合切向报告当成直接法向证书。

"""仅用标准库审计并汇总配对/母面验证器的球面和圆环面测量。"""
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
        if r["backend"] == "mother" and "joint_tangent_bound_deg" in r:
            assert float(r["tangent_u_bound_deg"]) < .1
            assert float(r["tangent_v_bound_deg"]) < .1
            assert float(r["joint_tangent_bound_deg"]) < .1
            assert float(r["joint_sine_squared_bound"]) < math.sin(math.radians(.1)) ** 2
            assert not r["joint_tangent_failure"]
        else:
            assert float(r["paired_tangent_bound"] if "paired_tangent_bound" in r else r["tangent_bound"]) < .1
        # 历史 CSV 保留其原有法向证书；当前母面验证记录没有直接法向约束，
        # 配对验证器后端仍保留其法向检查，不能混用两个合同。
        if "normal_bound" in r:
            assert float(r["normal_bound"]) < .1
        elif r["backend"] == "paired":
            assert float(r["paired_normal_bound"]) < .1
        else:
            assert not r["paired_normal_bound"]
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
    joint = "joint_tangent_bound_deg" in rows[0]
    summary = dict(schema=2 if joint else 1, input_rows=len(rows), runs_per_case_backend=3,
        contract=("mother distance with floating-point guard plus interval boundary position and joint U/V tangent budget"
                  if joint else "mother distance with floating-point guard plus interval boundary position/tangent"),
        joint_tangent_budget=joint,
        mother_normal_checked="normal_bound" in rows[0],
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
