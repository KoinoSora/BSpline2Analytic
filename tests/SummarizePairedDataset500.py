"""Audit all 500 new-only cases and aggregate exclusive validator modules.

The paired validator currently accepts polynomial input only. Rational cases
exercise Recognize separately; neither uncalled modules nor missing candidates
are averaged in as zero-cost successful validation.
"""
import argparse
import collections
import csv
import json
import math
from pathlib import Path
import statistics

MODULES = ("prepare_ms", "surface_ms", "subdivide_ms", "boundary_ms", "report_ms")
KINDS = ("Plane", "Cylinder", "Cone", "Sphere", "Torus")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def distribution(values):
    values = sorted(values)
    if not values:
        return {"count": 0, "total_ms": None, "mean_ms": None, "median_ms": None,
                "p95_ms": None, "max_ms": None, "under_10ms": 0}
    return {"count": len(values), "total_ms": sum(values),
            "mean_ms": statistics.mean(values), "median_ms": statistics.median(values),
            "p95_ms": values[math.ceil(.95 * len(values)) - 1], "max_ms": max(values),
            "under_10ms": sum(v < 10 for v in values)}


def paired_group(rows):
    measured = [r for r in rows if r["validation_ran"] == "1"]
    totals = {k: sum(float(r[k]) for r in measured) for k in MODULES}
    total = sum(totals.values())
    return {"cases": len(rows), "statuses": dict(collections.Counter(r["new_status"] for r in rows)),
            "measured": len(measured),
            "modules": {k: {"total_ms": value if measured else None,
                             "mean_ms": value / len(measured) if measured else None,
                             "percent": 100 * value / total if total else None}
                        for k, value in totals.items()},
            "validation": distribution([float(r["validation_total_ms"]) for r in measured]),
            # The historical driver does not finish this external stopwatch
            # when classification throws. Do not treat those zeros as timings.
            "classification_compute": distribution([float(r["classify_ms"]) for r in rows
                                                       if r["new_status"] != "ClassificationFailed"]),
            "classification_compute_unrecorded": sum(r["new_status"] == "ClassificationFailed" for r in rows),
            "classification_plus_first_validation": distribution(
                [float(r["classify_ms"]) + float(r["new_ms"]) for r in measured]),
            "successful_classification_plus_first_validation": distribution(
                [float(r["classify_ms"]) + float(r["new_ms"]) for r in measured
                 if r["new_status"] == "Passed"]),
            "surface_nodes": sum(int(r["surface_nodes"]) for r in measured),
            "boundary_nodes": sum(int(r["boundary_nodes"]) for r in measured)}


def exact_group(rows):
    return {"cases": len(rows), "statuses": dict(collections.Counter(r["exact_status"] for r in rows)),
            "recognition": distribution([float(r["exact_ms"]) for r in rows]),
            "successful_recognition": distribution([float(r["exact_ms"]) for r in rows
                                                      if r["exact_status"] == "Success"]),
            "classification_compute_total_ms": sum(float(r["exact_classify_compute_ms"]) for r in rows),
            "classification_validation_total_ms": sum(float(r["exact_classify_validate_ms"]) for r in rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.csv.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    require(len(rows) == 500 and {int(r["id"]) for r in rows} == set(range(500)),
            "Expected 500 unique cases, IDs 0..499")
    require(collections.Counter(r["type"] for r in rows) == {k: 100 for k in KINDS},
            "Surface-type quota mismatch")
    require(sum(int(r["rational"]) for r in rows) == 200, "Rational quota mismatch")
    require(sum(int(r["collapsed"]) != 0 for r in rows) == 40, "Collapsed quota mismatch")
    for r in rows:
        require(r["old_status"] == "NotRun" and float(r["old_ms"]) == 0,
                f"Legacy validator was invoked at {r['id']}")
        for k in MODULES + ("validation_total_ms", "classify_ms", "new_ms", "exact_ms"):
            require(math.isfinite(float(r[k])) and float(r[k]) >= 0, f"Invalid {k} at {r['id']}")
        require(math.isclose(sum(float(r[k]) for k in MODULES), float(r["validation_total_ms"]),
                             rel_tol=1e-8, abs_tol=1e-8), f"Timing mismatch at {r['id']}")
        if r["rational"] == "1":
            require(r["validation_ran"] == "0" and r["new_status"] == "Unsupported",
                    "Rational input counted as paired validation")
        if r["new_status"] == "Passed":
            require(r["validation_ran"] == "1" and r["collapsed"] == "0", "Invalid success")
            for k, limit in (("geometry_bound", 1e-5), ("tangent_deg", .1), ("normal_deg", .1)):
                require(math.isfinite(float(r[k])) and float(r[k]) < limit, f"Tolerance {k} violated")
    polynomial = [r for r in rows if r["rational"] == "0"]
    rational = [r for r in rows if r["rational"] == "1"]
    summary = {"source": str(args.csv), "cases": len(rows),
               "scope": "Polynomial first recovered candidate only; rational Recognize separately; no optimization or legacy projection validation",
               "module_denominator": "Cases with validation_ran=1, including unsuccessful attempts",
               "polynomial": paired_group(polynomial), "rational": exact_group(rational),
               "by_type": {}, "regular_failures": [], "degenerate_results": []}
    for kind in KINDS:
        selected = [r for r in rows if r["type"] == kind]
        summary["by_type"][kind] = {
            "cases": len(selected),
            "polynomial": paired_group([r for r in selected if r["rational"] == "0"]),
            "polynomial_passed": paired_group([r for r in selected if r["new_status"] == "Passed"]),
            "rational": exact_group([r for r in selected if r["rational"] == "1"])}
    for r in rows:
        status = r["exact_status"] if r["rational"] == "1" else r["new_status"]
        item = {k: r[k] for k in ("id", "type", "rational", "message")}
        item["status"] = status
        if int(r["collapsed"]):
            summary["degenerate_results"].append(item)
        elif status not in ("Passed", "Success"):
            summary["regular_failures"].append(item)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
                                               encoding="utf-8")
    with (args.output / "modules-by-type.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("type", "polynomial_cases", "measured_cases", "module",
                                                     "total_ms", "mean_ms", "percent"))
        writer.writeheader()
        for kind, group in summary["by_type"].items():
            g = group["polynomial"]
            for module, values in g["modules"].items():
                writer.writerow({"type": kind, "polynomial_cases": g["cases"], "measured_cases": g["measured"],
                                 "module": module, **values})
    print(json.dumps({"cases": len(rows), "polynomial": summary["polynomial"],
                      "rational": summary["rational"], "regular_failures": summary["regular_failures"]},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
