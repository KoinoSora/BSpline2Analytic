"""Run specialized Dataset500 cases with a per-case watchdog and no warm-up.

The watchdog includes loading/certification; measured stages come only from the
C++ timers. Missing measurements remain empty, never zero. Every selected case
is retained, including collapsed inputs and failures.
"""
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import time

MODULES = ("classify_compute_ms", "classify_validate_ms", "optimize_compute_ms", "optimize_validate_ms")
FIELDS = ("id", "pair", "type", "source_kind", "path", "controls", "degree_u", "degree_v",
          "collapsed_mask", "fixture_bound", "status", "classify_ms", "optimize_ms", *MODULES,
          "classify_optimize_ms", "total_ms", "geometry_bound", "tangent_deg", "normal_deg",
          "classify_compute_calls", "classify_validate_calls", "optimize_compute_calls",
          "optimize_validate_calls", "opt_iterations", "timing_schema", "fixture_qualified",
          "watchdog_elapsed_s", "message")


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("exe", "data", "sdk", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--timeout", type=float, default=20)
    parser.add_argument("--sphere-polynomial-only", action="store_true",
                        help="Select all 50 nonrational spheres, including 10 collapsed cases")
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be finite and positive")
    manifest = read(args.data / "manifest.csv")
    if [int(r["id"]) for r in manifest] != list(range(500)):
        raise ValueError("Expected the complete ordered 500-case manifest")
    if args.sphere_polynomial_only:
        manifest = [r for r in manifest if r["type"] == "Sphere" and r["source_kind"] == "refitted"]
        if [int(r["id"]) for r in manifest] != list(range(301, 400, 2)):
            raise ValueError("Expected the complete 50-case nonrational sphere subset")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    parts = args.output.with_suffix(".parts")
    if parts.exists():
        raise ValueError("Use a fresh output path; existing part results are never reused")
    parts.mkdir()
    env = {k.upper(): v for k, v in os.environ.items()}
    env["PATH"] = str(args.sdk.resolve() / "x64-win/bin") + os.pathsep + env.get("PATH", "")
    results = []
    with args.output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, extrasaction="ignore")
        writer.writeheader()
        for fixture in manifest:
            id_ = int(fixture["id"])
            part = parts / f"case-{id_:03}.csv"
            log = part.with_suffix(".log")
            start = time.perf_counter()
            code, timeout = None, False
            with log.open("wb") as capture:
                try:
                    code = subprocess.run([str(args.exe.resolve()), "--case", str(args.data.resolve()),
                                           str(part.resolve()), str(id_)], cwd=args.sdk.resolve(),
                                          env=env, stdout=capture, stderr=subprocess.STDOUT,
                                          timeout=args.timeout,
                                          creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0)).returncode
                except subprocess.TimeoutExpired:
                    timeout = True
            elapsed = time.perf_counter() - start
            rows = read(part) if part.exists() and part.stat().st_size else []
            began = f"BEGIN_CONVERT {id_}" in log.read_text(encoding="utf-8", errors="replace")
            if len(rows) == 1 and rows[0]["id"] == str(id_) and not timeout and code in (0, 1):
                row = rows[0]
                for key in MODULES:
                    if not math.isfinite(float(row[key])) or float(row[key]) < 0:
                        raise ValueError(f"Invalid timing at {id_}")
                if not math.isclose(sum(float(row[k]) for k in MODULES),
                                    float(row["classify_optimize_ms"]), abs_tol=1e-8, rel_tol=1e-8):
                    raise ValueError(f"Timing partition mismatch at {id_}")
            else:
                row = {k: fixture[k] for k in ("id", "pair", "type", "source_kind", "degree_u", "degree_v", "collapsed_mask")}
                row.update(controls=int(fixture["nu"]) * int(fixture["nv"]),
                           path="polynomial" if fixture["source_kind"] == "refitted" or fixture["type"] == "Plane" else "rational",
                           status="ProcessTimeout" if timeout else "ProcessFailure",
                           message=f"External watchdog: timeout={timeout}, child_exit={code}, conversion_started={began}")
            row.update(fixture_qualified=int(began), watchdog_elapsed_s=elapsed)
            expected = "Degenerate" if int(fixture["collapsed_mask"]) else "Success"
            if row["status"] == expected and code != 0:
                row.update(status="ProcessFailure", message=f"Result followed by child exit={code}")
            writer.writerow(row)
            stream.flush()
            results.append(row)
            print(f"{id_}: {row['status']} ({elapsed:.3f}s including loading/audits)", flush=True)
    unexpected = [int(r["id"]) for r in results
                  if r["status"] != ("Degenerate" if int(r["collapsed_mask"]) else "Success")]
    summary = {"cases": len(manifest), "warmup": False, "watchdog_seconds": args.timeout,
               "statuses": dict(collections.Counter(r["status"] for r in results)),
               "unexpected_ids": unexpected,
               "qualified_inputs": sum(int(r["fixture_qualified"]) for r in results),
               "successful_regular_under_10ms": sum(r["status"] == "Success" and
                   float(r["classify_optimize_ms"]) < 10 for r in results),
               "regular_denominator": sum(not int(r["collapsed_mask"]) for r in manifest)}
    args.output.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return bool(unexpected)


if __name__ == "__main__":
    raise SystemExit(main())
