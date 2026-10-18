# 文件职责：历史一般 ConvertGeneral 记录的统计工具。
# 检查清单身份和候选计数；旧数据结果不能直接代表当前 V2 合同。

"""仅用标准库检查成员集合并汇总历史 ConvertGeneral 运行记录。"""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+", type=Path)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--extras", action="store_true")
    selection.add_argument("--positives", action="store_true")
    parser.add_argument("--merge", type=Path)
    parser.add_argument("--configuration", choices=("64-full", "8-full", "8-staged"),
                        help="Select one configuration from --compare-positive output")
    args = parser.parse_args()
    rows = []
    for path in args.csv:
        with path.open(encoding="utf-8-sig", newline="") as source:
            rows.extend(csv.DictReader(source))
    if args.configuration:
        rows = [r for r in rows if r.get("configuration") == args.configuration]
    rows.sort(key=lambda r: int(r["id"]))
    count = 60 if args.extras else 500 if args.positives else 600
    if [int(r["id"]) for r in rows] != list(range(count)):
        raise SystemExit(f"Expected exactly IDs 0..{count-1}, without gaps or duplicates")
    if any(int(r["extra"]) != args.extras or r["fixture_certified"] != "1" for r in rows):
        raise SystemExit("Dataset identity or independent input certificate mismatch")
    if not args.extras:
        expected = ["feasible"] * 500 + ([] if args.positives else ["infeasible"] * 100)
        if [r["group"] for r in rows] != expected:
            raise SystemExit("Fixture groups do not match the requested dataset")
        if sum(int(r["rational"]) for r in rows[:500]) != 250:
            raise SystemExit("Positive rational quota mismatch")
        if not args.positives and sum(int(r["rational"]) for r in rows[500:]) != 50:
            raise SystemExit("Negative rational quota mismatch")
        if any(not 50 <= int(r["controls"]) <= 200 for r in rows[:500]):
            raise SystemExit("Positive control budget violated")
    if args.merge:
        with args.merge.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(destination, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
    print(f"{count} complete rows, all input certificates verified")
    if all("score_evaluations" in row for row in rows):
        if any(int(r["raw_candidates"]) != int(r["duplicate_candidates"]) +
               int(r["score_evaluations"]) for r in rows):
            raise SystemExit("Candidate deduplication/scoring counters are inconsistent")
        raw = sum(int(r["raw_candidates"]) for r in rows)
        duplicates = sum(int(r["duplicate_candidates"]) for r in rows)
        scored = sum(int(r["score_evaluations"]) for r in rows)
        print(f"Generated candidates: {raw}; duplicates skipped: {duplicates}; score calls: {scored}")
    if all("full_scores" in row for row in rows):
        if any(not 0 <= int(r["completed_scores"]) <= int(r["full_scores"]) <=
               int(r["score_evaluations"]) for r in rows):
            raise SystemExit("Progressive scoring counters are inconsistent")
        for field in ("seed_attempts", "coarse_scores", "full_scores", "completed_scores",
                      "deferred_candidates", "backup_candidates", "interior_point_evaluations",
                      "boundary_point_evaluations"):
            print(f"{field}: {sum(int(r[field]) for r in rows)}")
    if all("torus_groups" in row for row in rows):
        if any(not 0 <= int(r["torus_candidates"]) <= int(r["torus_axes"]) <=
               2 * int(r["torus_groups"]) or
               int(r["torus_candidates"]) + int(r["torus_circle_failures"]) +
               int(r["torus_invalid_radii"]) != int(r["torus_axes"]) or
               2 * int(r["torus_dual_groups"]) > int(r["torus_candidates"]) for r in rows):
            raise SystemExit("Direct torus construction counters are inconsistent")
        for field in ("torus_groups", "torus_rank_deficient", "torus_axis_degenerate",
                      "torus_no_real_axis", "torus_axes", "torus_rejected_axes",
                      "torus_circle_failures", "torus_invalid_radii", "torus_candidates",
                      "torus_dual_groups"):
            print(f"{field}: {sum(int(r[field]) for r in rows)}")
    print("Stage totals (seconds): " + ", ".join(
        f"{field}={sum(float(r[field]) for r in rows)/1000:.6f}"
        for field in ("classify_ms", "optimize_ms", "validate_ms", "total_ms")))
    groups = defaultdict(list)
    for row in rows:
        groups[(row["group"], row["source_type"])].append(row)
    print("| Group | Type | Cases | Success | Expected checks passed | Total mean ms | Total P95 ms | Max classify+optimize ms |")
    print("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |")
    for (group, kind), values in sorted(groups.items()):
        times = sorted(float(v["total_ms"]) for v in values)
        solve = max(float(v["classify_ms"]) + float(v["optimize_ms"]) for v in values)
        print(f"| {group} | {kind} | {len(values)} | {sum(int(v['success']) for v in values)} | "
              f"{sum(int(v['passed']) for v in values)} | {sum(times)/len(times):.3f} | "
              f"{times[math.ceil(.95*len(times))-1]:.3f} | {solve:.3f} |")
    failures = defaultdict(list)
    for row in rows:
        if row["passed"] != "1":
            reason = row["pcurve_failure"]
            if row["success"] == "1" and row["output_audit"] != "1":
                reason = "IndependentOutputAudit"
            failures[(row["group"], reason)].append(row["file"])
    print("\nUnmet expectations:")
    for (group, reason), files in sorted(failures.items()):
        print(f"- {group} / {reason}: {len(files)}; first={files[0]}")
    print(f"Checks passed: {sum(int(r['passed']) for r in rows)}/{count}")
    false_accepts = [r["file"] for r in rows if r["group"] == "infeasible" and r["success"] == "1"]
    print("Certified negatives: not run" if args.positives else f"Certified-negative false accepts: {len(false_accepts)}")
    for filename in false_accepts:
        print(filename)
    if failures:
        raise SystemExit("Unmet test expectations")


if __name__ == "__main__":
    main()
