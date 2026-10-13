"""Validate the complete extreme batch and write reproducible grouped timings.

Exit zero means the CSVs are structurally valid, not that all conversions pass.
No SDK, conversion, or old dataset is invoked by this script.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean


def read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def finite_max(rows, key):
    values = [float(r[key]) for r in rows]
    return max((x for x in values if math.isfinite(x)), default=None)


def summarize(manifest, results, modules, output):
    source, rows, probes = read(manifest), read(results), read(modules)
    expected = list(range(200))
    for label, data in (("manifest", source), ("results", rows)):
        if [int(row["id"]) for row in data] != expected:
            raise ValueError(f"{label}: expected exactly one ordered row per ID 0..199")
    if len(probes) != 7:
        raise ValueError("Expected five biased-seed probes and two angle-rejection probes")
    for fixture, row in zip(source, rows):
        for key in ("id", "group", "type", "source_kind", "controls", "rational"):
            if fixture[key] != row[key]:
                raise ValueError(f"Fixture/result mismatch: {row['id']} {key}")
        times = [float(row[key]) for key in
                 ("classify_ms", "optimize_ms", "validate_ms", "total_ms", "classify_optimize_ms")]
        if row["status"] in ("process_timeout", "process_failure"):
            if not all(math.isnan(x) for x in times):
                raise ValueError("Incomplete process must leave stage timings unknown")
            continue
        if not all(math.isfinite(x) and x >= 0 for x in times):
            raise ValueError(f"Invalid timing: {row['id']}")
        if not math.isclose(times[0] + times[1], times[4], rel_tol=1e-12, abs_tol=1e-9):
            raise ValueError(f"Combined timing is not a per-case sum: {row['id']}")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(list)
    for row in rows:
        groups["all", "all"].append(row)
        groups["scenario", row["group"]].append(row)
        groups["type_source", row["type"] + "/" + row["source_kind"]].append(row)
    fields = ["dimension", "group", "count", "same_type", "alternate_type", "failed",
              "sampled_violation", "qualified", "sampled_only_angles", "unknown_timings",
              "classify_mean_ms", "optimize_mean_ms", "validate_mean_ms", "total_mean_ms",
              "classify_optimize_mean_ms", "classify_optimize_max_ms", "slowest_id"]
    aggregates = []
    for (dimension, name), members in sorted(groups.items()):
        counts = Counter(r["status"] for r in members)
        timed = [r for r in members if math.isfinite(float(r["total_ms"]))]
        slowest = max(timed, key=lambda r: float(r["classify_optimize_ms"]), default=None)
        item = dict(dimension=dimension, group=name, count=len(members),
                    same_type=counts["same_type_success"], alternate_type=counts["alternate_type_success"],
                    failed=len(members)-counts["same_type_success"]-counts["alternate_type_success"],
                    sampled_violation=counts["success_with_sampled_violation"],
                    qualified=sum(r["qualified"] == "1" for r in members),
                    sampled_only_angles=sum(r["angle_evidence"] == "sampled_only" for r in members),
                    unknown_timings=len(members)-len(timed),
                    classify_optimize_max_ms=float(slowest["classify_optimize_ms"]) if slowest else None,
                    slowest_id=slowest["id"] if slowest else None)
        for stage in ("classify", "optimize", "validate", "total", "classify_optimize"):
            item[stage+"_mean_ms"] = mean(float(r[stage+"_ms"]) for r in timed) if timed else None
        aggregates.append(item)
    with (output / "groups.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(aggregates)
    summary = {
        "counts": dict(Counter(r["status"] for r in rows)),
        "modules_passed": sum(r["passed"] == "1" for r in probes),
        "source_distance_sampled_max": finite_max(rows, "source_distance_sampled"),
        "source_tangent_sampled_deg_max": finite_max(rows, "source_tangent_sampled_deg"),
        "source_normal_sampled_deg_max": finite_max(rows, "source_normal_sampled_deg"),
        "pair_bound_max": max(float(r["pair_bound"]) for r in rows),
        "slowest_classify_optimize": sorted((r for r in rows if math.isfinite(float(r["total_ms"]))), key=lambda r: float(r["classify_optimize_ms"]), reverse=True)[:10],
        "slowest_total": sorted((r for r in rows if math.isfinite(float(r["total_ms"]))), key=lambda r: float(r["total_ms"]), reverse=True)[:10],
        "input_hashes": {Path(p).name: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                         for p in (manifest, results, modules)},
        "groups": aggregates,
    }
    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, list)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest")
    parser.add_argument("results")
    parser.add_argument("modules")
    parser.add_argument("output_directory")
    args = parser.parse_args()
    summarize(args.manifest, args.results, args.modules, args.output_directory)
