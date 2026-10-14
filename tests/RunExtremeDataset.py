"""Run only Extreme200, isolating each SDK call in a bounded child process.

Timeout/crash rows contain NaN (unknown) stage timings, never invented zeros.
--resume preserves completed rows from this same batch; CTest starts a fresh run.
"""
import argparse
import csv
import math
import subprocess
import time
from pathlib import Path

FIELDS = ("id,group,type,source_kind,qualified,controls,rational,pair_bound,"
          "source_distance_sampled,source_tangent_sampled_deg,source_normal_sampled_deg,"
          "status,path,output_type,classify_ms,optimize_ms,classify_optimize_ms,validate_ms,"
          "total_ms,geometry_bound,boundary_bound,output_distance_sampled,"
          "output_tangent_sampled_deg,output_normal_sampled_deg,angle_evidence,"
          "iterations,candidates,pieces,attempts,message,watchdog_elapsed_s").split(",")


def read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--last", type=int, default=200, help="exclusive end ID")
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be finite and positive")
    if not 0 <= args.first < args.last <= 200:
        parser.error("Expected 0 <= first < last <= 200")
    manifest = read(args.data / "manifest.csv")
    if [int(r["id"]) for r in manifest] != list(range(200)):
        raise ValueError("Expected exactly 200 ordered fixtures")
    manifest = manifest[args.first:args.last]
    existing = read(args.output) if args.resume and args.output.exists() else []
    if len(existing) > len(manifest) or [int(r["id"]) for r in existing] != list(range(args.first, args.first + len(existing))):
        raise ValueError("Resume requires a contiguous prefix of the selected range")
    for row, fixture in zip(existing, manifest):
        if any(row[k] != fixture[k] for k in ("id", "group", "type", "source_kind")):
            raise ValueError("Resume prefix belongs to a different matrix")
    partdir = args.output.with_suffix(".parts")
    partdir.mkdir(parents=True, exist_ok=True)
    failed = sum(r["status"] not in ("same_type_success", "alternate_type_success") for r in existing)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS)
        writer.writeheader()
        writer.writerows(existing)
        stream.flush()
        for fixture in manifest[len(existing):]:
            index = int(fixture["id"])
            part = partdir / f"case-{index:03}.csv"
            # Remove only this runner's old per-case output, never input geometry.
            part.unlink(missing_ok=True)
            log = part.with_suffix(".log")
            start = time.perf_counter()
            timeout = False
            with log.open("wb") as capture:
                try:
                    run = subprocess.run([str(args.exe.resolve()), "--run", str(args.data.resolve()),
                                          str(part.resolve()), str(index), str(index+1)],
                                         stdout=capture, stderr=subprocess.STDOUT, timeout=args.timeout,
                                         creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                    code = run.returncode
                except subprocess.TimeoutExpired:
                    timeout = True
                    code = None
            elapsed = time.perf_counter() - start
            completed = read(part) if part.exists() and part.stat().st_size else []
            if len(completed) == 1 and completed[0]["id"] == str(index) and not timeout:
                row = completed[0]
                if code not in (0, 1):
                    row["status"] = "process_failure_after_result"
                    row["message"] += f"; child exit={code}"
            else:
                row = {key: "nan" for key in FIELDS}
                for key in ("id", "group", "type", "source_kind", "controls", "rational", "pair_bound"):
                    row[key] = fixture[key]
                began = f"BEGIN_CONVERT {index}" in log.read_text(encoding="utf-8", errors="replace")
                row.update(qualified="1" if began else "0", status="process_timeout" if timeout else "process_failure",
                           path="unknown", output_type="Unknown", angle_evidence="not_checked",
                           attempts="incomplete", message=f"External watchdog; child exit={code}; conversion_started={began}")
            row["watchdog_elapsed_s"] = elapsed
            writer.writerow(row)
            stream.flush()
            failed += row["status"] not in ("same_type_success", "alternate_type_success")
            print(f"{index}: {row['status']} ({elapsed:.3f}s including loading and audits)", flush=True)
    print(f"Completed {len(manifest)} records [{args.first},{args.last}), failed/unresolved={failed}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
