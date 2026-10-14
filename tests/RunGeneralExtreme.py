"""Run only the 320 general extreme inputs, one bounded SDK process per input."""
import argparse
import csv
import hashlib
import math
import os
import subprocess
import time
from pathlib import Path
from SummarizeGeneralExtreme import audit_data

FIELDS = ("id,group,variant,type,rational,qualified,source_bound,status,output_type,output_bound,"
          "conversion_status,pcurve_failure,classify_ms,optimize_ms,validate_ms,total_ms,geometry_bound,"
          "boundary_bound,attempts,iterations,validation_pieces,chart_attempts,chart_iterations,seed_groups,"
          "raw_candidates,duplicate_candidates,full_scores,torus_groups,torus_rank_deficient,"
          "torus_no_real_axis,torus_candidates,torus_dual_groups,message,input_sha256,watchdog_seconds").split(",")
SUCCESS = {"same_type_success", "alternate_type_success"}


def read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--sdk", required=True, type=Path)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--configuration", choices=("Debug", "Release", "RelWithDebInfo", "MinSizeRel"), default="Release")
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--last", type=int, default=320)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.first < args.last <= 320 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("Invalid range or timeout")
    manifest, _ = audit_data(args.data)
    if [int(r["id"]) for r in manifest] != list(range(320)):
        raise ValueError("Expected complete GeneralExtreme320 manifest")
    selected = manifest[args.first:args.last]
    existing = read(args.output) if args.resume and args.output.exists() else []
    if len(existing) > len(selected) or [int(r["id"]) for r in existing] != list(range(args.first, args.first+len(existing))):
        raise ValueError("Resume needs a contiguous prefix of the selected range")
    for row, fixture in zip(existing, selected):
        digest = hashlib.sha256((args.data / fixture["file"]).read_bytes()).hexdigest()
        if row["input_sha256"] != digest:
            raise ValueError("Input changed since the previous run")
    env = dict(os.environ)
    suffix = {"Debug": "d", "RelWithDebInfo": "w"}.get(args.configuration, "")
    env["PATH"] = str(args.sdk.resolve() / ("x64-win/bin"+suffix)) + ";" + env.get("PATH", "")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    parts = args.output.with_suffix(".parts")
    parts.mkdir(exist_ok=True)
    failures = sum(r["status"] not in SUCCESS for r in existing)
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS)
        writer.writeheader()
        writer.writerows(existing)
        stream.flush()
        for fixture in selected[len(existing):]:
            index = int(fixture["id"])
            output = parts / f"case-{index:03}.csv"
            output.unlink(missing_ok=True)  # Only this runner's own output.
            log = output.with_suffix(".log")
            start = time.perf_counter()
            timeout = False
            with log.open("wb") as capture:
                try:
                    run = subprocess.run([str(args.exe.resolve()), "--run", str(args.data.resolve()),
                                          str(output.resolve()), str(index)], cwd=args.sdk.resolve(), env=env,
                                         stdout=capture, stderr=subprocess.STDOUT, timeout=args.timeout,
                                         creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
                    code = run.returncode
                except subprocess.TimeoutExpired:
                    timeout, code = True, None
            elapsed = time.perf_counter()-start
            rows = read(output) if output.exists() and output.stat().st_size else []
            content = log.read_text(encoding="utf-8", errors="replace")
            phase = "input" if f"BEGIN_CONVERT {index}" not in content else "convert" if f"END_CONVERT {index}" not in content else "audit"
            if len(rows) == 1 and rows[0]["id"] == str(index):
                row = rows[0]
                if timeout or code not in (0, 1):
                    row["status"] = phase + ("_timeout" if timeout else "_process_failure")
                    row["message"] += f"; child exit={code}"
            else:
                row = {key: "nan" for key in FIELDS}
                for key in ("id", "group", "variant", "type", "rational"):
                    row[key] = fixture[key]
                row.update(qualified="1" if phase != "input" else "0", source_bound=fixture["source_bound"],
                           status=phase+("_timeout" if timeout else "_process_failure"), output_type="Unknown",
                           pcurve_failure="not_available", message=f"phase={phase}; child exit={code}")
            row["input_sha256"] = hashlib.sha256((args.data/fixture["file"]).read_bytes()).hexdigest()
            row["watchdog_seconds"] = elapsed
            writer.writerow(row)
            stream.flush()
            failures += row["status"] not in SUCCESS
            print(f"{index}: {row['status']} ({elapsed:.3f}s incl. certificates)", flush=True)
    print(f"Completed {len(selected)}; failed/unresolved={failures}", flush=True)
    return int(failures != 0)


if __name__ == "__main__":
    raise SystemExit(main())
