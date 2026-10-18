# 文件职责：历史专用转换阶段计时的汇总工具。
# 校验细分阶段与总计关系，不对已废弃运行入口进行调用。

"""校验并汇总 specialized-v2 互斥阶段计时划分。"""
import argparse
import collections
import csv
import json
import math
from pathlib import Path

MODULES = ("classify_compute_ms", "classify_validate_ms",
           "optimize_compute_ms", "optimize_validate_ms")
METRICS = MODULES + ("classify_ms", "optimize_ms", "classify_optimize_ms", "other_ms", "total_ms")


def require(ok, message):
    if not ok:
        raise ValueError(message)


def summarize(rows):
    result = {"count": len(rows), "statuses": dict(collections.Counter(r["status"] for r in rows))}
    for name in METRICS:
        ordered = sorted(rows, key=lambda r: (r[name], r["id"]))
        values = [r[name] for r in ordered]
        result[name] = {"sum": math.fsum(values), "mean": math.fsum(values) / len(values),
                        "p95": values[math.ceil(.95 * len(values)) - 1], "max": values[-1],
                        "max_id": min(r["id"] for r in rows if r[name] == values[-1])}
    result["iterations"] = sum(r["opt_iterations"] for r in rows)
    result["optimized_cases"] = sum(r["optimize_compute_calls"] > 0 for r in rows)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.csv.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 500, "Expected 500 rows")
    integers = ("id", "collapsed_mask", "controls", "opt_iterations", "classify_compute_calls",
                "classify_validate_calls", "optimize_compute_calls", "optimize_validate_calls")
    for r in rows:
        require(r["timing_schema"] == "specialized-v2", "Wrong timing schema")
        for name in integers:
            r[name] = int(r[name])
        for name in METRICS + ("validate_ms",):
            r[name] = float(r[name])
            require(math.isfinite(r[name]) and r[name] >= -1e-8, f"Invalid {name} at {r['id']}")
        require(r["path"] in ("rational", "polynomial"), "Non-specialized path")
        expected = "Degenerate" if r["collapsed_mask"] else "Success"
        r["expected_status"] = expected
        for phase, compute, validate in (("classify_ms", MODULES[0], MODULES[1]),
                                         ("optimize_ms", MODULES[2], MODULES[3])):
            require(math.isclose(r[phase], r[compute] + r[validate], rel_tol=1e-8, abs_tol=1e-8),
                    f"Phase partition mismatch at {r['id']}")
        require(r["validate_ms"] == 0, "Validation counted outside its owning phase")
        require(math.isclose(r["classify_optimize_ms"], sum(r[k] for k in MODULES),
                             rel_tol=1e-8, abs_tol=1e-8), "Evaluation timing mismatch")
        require(math.isclose(r["total_ms"], r["classify_optimize_ms"] + r["other_ms"],
                             rel_tol=1e-8, abs_tol=1e-8), "Total partition mismatch")
        for module in MODULES:
            count = r[module.removesuffix("_ms") + "_calls"]
            require(count or r[module] == 0, "Unused module has nonzero time")
        if r["status"] == "Success":
            require(r["classify_compute_calls"] > 0 and r["classify_validate_calls"] > 0,
                    "Successful specialized result has no classification/validation calls")
        if r["path"] == "rational":
            require(r["optimize_ms"] == 0 and r["opt_iterations"] == 0 and
                    r["optimize_compute_calls"] == 0 and r["optimize_validate_calls"] == 0,
                    "Exact path optimized")
    require(sorted(r["id"] for r in rows) == list(range(500)), "Missing or duplicate IDs")
    require(collections.Counter(r["path"] for r in rows) == {"rational": 200, "polynomial": 300},
            "Wrong representation distribution")
    failures = [r["id"] for r in rows if r["status"] != r["expected_status"]]
    groups = collections.defaultdict(list)
    for r in rows:
        kind = "degenerate" if r["collapsed_mask"] else "regular"
        for group in (kind, f"{kind}/{r['path']}", f"{kind}/{r['type']}/{r['path']}",
                      f"{kind}/bucket-{(r['id'] % 100) // 20}/{r['path']}"):
            groups[group].append(r)
    summary = {"timing_schema": "specialized-v2", "input": args.csv.name,
               "all": summarize(rows), "unexpected_ids": failures,
               "groups": {k: summarize(v) for k, v in sorted(groups.items())}}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    flat = []
    for name, group in summary["groups"].items():
        item = {"group": name, "count": group["count"], "iterations": group["iterations"],
                "optimized_cases": group["optimized_cases"]}
        for metric in METRICS:
            item.update({f"{metric}_{key}": value for key, value in group[metric].items()})
        flat.append(item)
    with (args.output / "groups.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    print(json.dumps({"count": len(rows), "statuses": summary["all"]["statuses"],
                      "unexpected_ids": failures, "totals_ms": {k: summary['all'][k]['sum'] for k in MODULES}}, indent=2))
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
