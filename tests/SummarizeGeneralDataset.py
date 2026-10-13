"""Validate membership and summarize direct ConvertGeneral runs (stdlib only)."""
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
    args = parser.parse_args()
    rows = []
    for path in args.csv:
        with path.open(encoding="utf-8-sig", newline="") as source:
            rows.extend(csv.DictReader(source))
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


if __name__ == "__main__":
    main()
