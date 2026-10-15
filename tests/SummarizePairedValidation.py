"""Audit standalone paired-validator CSVs and compare the unchanged production path."""
import argparse
import collections
import csv
import json
import math
from pathlib import Path

PREPARATION_DETAILS = ("extraction_ms", "source_cache_ms", "local_transform_ms")


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def aggregate(data):
    selected = [r for r in data if r["old_status"] == r["new_status"] == "Passed"]
    old = sum(float(r["old_ms"]) for r in selected)
    new = sum(float(r["new_ms"]) for r in selected)
    stages = {k: sum(float(r[k]) for r in selected) for k in
              ("prepare_ms", "surface_ms", "subdivide_ms", "boundary_ms", "report_ms")}
    return {"rows": len(data), "statuses": dict(collections.Counter(r["new_status"] for r in data)),
            "matched_passes": len(selected), "old_ms": old, "new_ms": new,
            "speedup": old / new if new else None,
            "reduction_percent": 100 * (1 - new / old) if old else None,
            "old_mean_ms": old / len(selected) if selected else None,
            "new_mean_ms": new / len(selected) if selected else None,
            "stages_ms": stages,
            "preparation_details_ms": {k: sum(float(r[k]) for r in selected)
                                       for k in PREPARATION_DETAILS
                                       if selected and k in selected[0]},
            "stages_percent": {k: 100 * v / sum(stages.values()) for k, v in stages.items()} if selected else {},
            "mapping_attempts": dict(collections.Counter(r["attempts"] for r in selected)),
            "max_geometry_bound": max((float(r["geometry_bound"]) for r in selected), default=None),
            "max_tangent_deg": max((float(r["tangent_deg"]) for r in selected), default=None),
            "max_normal_deg": max((float(r["normal_deg"]) for r in selected), default=None)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--baseline", type=Path)
    args = parser.parse_args()
    root = args.directory
    runs = [rows(root / f"sphere-run{i}.csv") for i in range(1, 4)]
    expected = set(range(301, 400, 2))
    for run in runs:
        require(len(run) == 50 and {int(r["id"]) for r in run} == expected, "Sphere case set mismatch")
        for r in run:
            require(r["type"] == "Sphere" and r["rational"] == "0", "Wrong sphere input kind")
            if r["new_status"] == "Passed":
                require(r["collapsed"] == "0", "Collapsed edge incorrectly accepted")
                for k, limit in (("geometry_bound", 1e-5), ("tangent_deg", .1), ("normal_deg", .1)):
                    require(math.isfinite(float(r[k])) and float(r[k]) < limit, "Tolerance violation")
    stable = ("old_status", "new_status", "geometry_bound", "tangent_deg", "normal_deg",
              "attempts", "u_segments", "v_segments", "common_cells", "surface_nodes", "boundary_nodes")
    stable += tuple(k for k in ("source_extractions", "cache_hits", "u_matrices", "v_matrices")
                    if k in runs[0][0])
    for left, right in zip(runs, runs[1:]):
        require(all(a[k] == b[k] for a, b in zip(left, right) for k in stable), "Unstable repeated certificate")
    data = rows(root / "dataset500.csv")
    require(len(data) == 500 and {int(r["id"]) for r in data} == set(range(500)), "Dataset case set mismatch")
    for r in [r for run in runs for r in run] + data:
        if all(k in r for k in PREPARATION_DETAILS):
            nested = sum(float(r[k]) for k in PREPARATION_DETAILS + ("nurbs_ms", "mapping_ms"))
            require(nested <= float(r["prepare_ms"]) + 1e-5, "Preparation sub-timers overlap")
    summary = {"scope": "First recovered candidate; no optimizer or production fallback is invoked",
               "sphere_runs": [aggregate(r) for r in runs],
               "sphere_all_runs": aggregate([r for run in runs for r in run]),
               "dataset500": aggregate(data),
               "dataset_by_type": {kind: aggregate([r for r in data if r["type"] == kind])
                                   for kind in sorted({r["type"] for r in data})}}
    means = []
    for trio in zip(*runs):
        mean = dict(trio[0])
        for key in ("old_ms", "new_ms", "classify_ms", "prepare_ms", "surface_ms", "subdivide_ms",
                    "boundary_ms", "report_ms", "nurbs_ms", "mapping_ms") + tuple(
                        k for k in PREPARATION_DETAILS if k in mean):
            mean[key] = sum(float(r[key]) for r in trio) / len(trio)
        mean["speedup"] = mean["old_ms"] / mean["new_ms"] if mean["new_ms"] else 0
        mean["old_classification_total_ms"] = mean["classify_ms"] + mean["old_ms"]
        mean["new_classification_total_ms"] = mean["classify_ms"] + mean["new_ms"]
        means.append(mean)
    with (root / "sphere-case-means.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(means[0]))
        writer.writeheader()
        writer.writerows(means)
    if args.baseline:
        before = {r["id"]: r for r in rows(args.baseline)}
        after = rows(root / "legacy-dataset500.csv")
        require(len(before) == len(after) == 500, "Incomplete production regression")
        keys = ("status", "geometry_bound", "tangent_deg", "normal_deg", "opt_iterations",
                "class_raw_candidates", "class_recoveries", "validation_pieces")
        differences = [{"id": r["id"], "field": k, "before": before[r["id"]][k], "after": r[k]}
                       for r in after for k in keys if before[r["id"]][k] != r[k]]
        summary["production_regression"] = {
            "differences": differences,
            "statuses": dict(collections.Counter(r["status"] for r in after))}
    (root / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    require(not summary.get("production_regression", {}).get("differences"), "Production behavior changed")


if __name__ == "__main__":
    main()
