"""Summarize all 500 fixtures, cones, or polynomial spheres; Python standard library only."""
import argparse
import csv
import math
import statistics
from collections import Counter, defaultdict


def classify_optimize_ms(row):
    """Sum stages for one case before taking a group maximum."""
    return float(row["classify_ms"]) + float(row["optimize_ms"])


def summarize(path, cone_only=False, sphere_polynomial_only=False):
    with open(path, newline="", encoding="utf-8-sig") as source:
        rows = list(csv.DictReader(source))
    if cone_only and sphere_polynomial_only:
        raise ValueError("Choose only one subset")
    expected = (set(range(301, 400, 2)) if sphere_polynomial_only else
                set(range(200, 300)) if cone_only else set(range(500)))
    if len(rows) != len(expected) or {int(r["id"]) for r in rows} != expected:
        raise ValueError("Expected exactly one result for every case " +
                         ("301,303,..399" if sphere_polynomial_only else
                          "200..299" if cone_only else "0..499"))
    if cone_only and any(r["type"] != "Cone" for r in rows):
        raise ValueError("Cone-only CSV contains a different surface type")
    if sphere_polynomial_only and any(r["type"] != "Sphere" or
                                     r["path"] != "polynomial" or
                                     r["source_kind"] != "refitted" for r in rows):
        raise ValueError("Polynomial-sphere CSV contains a different type or input kind")
    regular = [r for r in rows if int(r["collapsed_mask"]) == 0]
    degenerate = [r for r in rows if int(r["collapsed_mask"]) != 0]
    print(f"Cases: {len(rows)}; regular: {len(regular)}; collapsed: {len(degenerate)}")
    print("Regular status:", dict(Counter(r["status"] for r in regular)))
    print("Collapsed status:", dict(Counter(r["status"] for r in degenerate)))
    print("Paths:", dict(Counter(r["path"] for r in rows)))
    refitted = [r for r in rows if r["source_kind"] == "refitted"]
    print("Max refit certificate:", max(float(r["fixture_bound"]) for r in refitted))
    print("Control count range:", min(int(r["controls"]) for r in rows),
          max(int(r["controls"]) for r in rows))

    def table(grouping, title):
        groups = defaultdict(list)
        for row in regular:
            groups[grouping(row)].append(row)
        print("\n" + title)
        print("| Group | N | Classification ms | Optimization ms | Validation ms | Total ms | Total P95 ms | Classification + optimization max ms | Max case |")
        print("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |")
        for key, group in sorted(groups.items()):
            means = [statistics.mean(float(r[k]) for r in group)
                     for k in ("classify_ms", "optimize_ms", "validate_ms", "total_ms")]
            values = sorted(float(r["total_ms"]) for r in group)
            p95 = values[math.ceil(.95 * len(values)) - 1]
            # Resolve equal maxima by case id so CSV row order cannot change the report.
            slowest = max(group, key=lambda r: (classify_optimize_ms(r), -int(r["id"])))
            print(f"| {key} | {len(group)} | " + " | ".join(
                f"{v:.4f}" for v in [*means, p95, classify_optimize_ms(slowest)])
                + f" | case-{int(slowest['id']):03}.b2a |")

    table(lambda r: r["type"] + "/" + r["path"], "By type and executed path (regular cases only)")
    limits = [(50, 79), (80, 109), (110, 139), (140, 169), (170, 200)]

    def band(row):
        count = int(row["controls"])
        lo, hi = next((lo, hi) for lo, hi in limits if lo <= count <= hi)
        return f"{lo:03}-{hi:03}/" + row["path"]

    table(band, "By control count and executed path (regular cases only)")
    for path_name in ("rational", "polynomial"):
        group = [r for r in degenerate if r["path"] == path_name]
        if not group:
            continue
        slowest = max(group, key=lambda r: (classify_optimize_ms(r), -int(r["id"])))
        print(f"Collapsed/{path_name}: n={len(group)}, mean total="
              f"{statistics.mean(float(r['total_ms']) for r in group):.4f} ms, "
              f"classification + optimization max={classify_optimize_ms(slowest):.4f} ms "
              f"(case-{int(slowest['id']):03}.b2a)")
    successful = [r for r in regular if r["status"] == "Success" and r["path"] == "polynomial"]
    for name in ("geometry_bound", "tangent_deg", "normal_deg"):
        if successful:
            print("Polynomial success max " + name + ":", max(float(r[name]) for r in successful))
    for row in sorted(regular, key=lambda r: float(r["total_ms"]), reverse=True)[:5]:
        print("Slow case:", row["id"], row["type"], row["path"], row["status"], row["total_ms"], "ms")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--cone-only", action="store_true")
    selection.add_argument("--sphere-polynomial-only", action="store_true")
    args = parser.parse_args()
    summarize(args.csv_path, args.cone_only, args.sphere_polynomial_only)
