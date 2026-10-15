"""Audit Dataset500 representation migration and record reproducible SHA-256 evidence.

This checks identity/representation preservation, not geometry. Run the C++ --audit
command as well for independent continuous source and approximation certificates.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    raw = path.read_bytes()
    t = raw.decode("utf-8").split()
    require(t[0] == "B2A_CASE_V1", f"Bad header: {path}")

    def surface(start):
        p, q, nu, nv, weights = map(int, t[start:start + 5])
        i = start + 5
        for _ in range(2):
            i += 1 + 2 * int(t[i])
        i += nu * nv * (3 + weights)
        return {"p": p, "q": q, "nu": nu, "nv": nv, "weights": weights,
                "tokens": t[start:i]}, i

    source, i = surface(26)
    target, i = surface(i)
    require(i == len(t), f"Truncated or trailing data: {path}")
    return {"metadata": t[1:9], "analytic": t[10:26], "source": source,
            "input": target, "sha256": hashlib.sha256(raw).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("original", type=Path)
    parser.add_argument("migrated", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows, grids, changed = [], [], []
    degrees = collections.Counter()
    with (args.migrated / "manifest.csv").open(encoding="utf-8-sig", newline="") as stream:
        manifest = list(csv.DictReader(stream))
    require([int(r["id"]) for r in manifest] == list(range(500)), "Manifest identity mismatch")
    for directory in (args.original, args.migrated):
        require(len(list(directory.glob("case-*.b2a"))) == 500, "Expected exactly 500 fixtures")
    for id_ in range(500):
        name = f"case-{id_:03}.b2a"
        old, new = read(args.original / name), read(args.migrated / name)
        require(old["metadata"] == new["metadata"] and int(new["metadata"][0]) == id_,
                f"Metadata changed: {name}")
        require(old["analytic"] == new["analytic"], f"Analytic truth/domain changed: {name}")
        require(old["source"] == new["source"], f"Saved source changed: {name}")
        a, b = old["input"], new["input"]
        require(all(int(manifest[id_][field]) == b[key] for field, key in
                    (("degree_u", "p"), ("degree_v", "q"), ("nu", "nu"), ("nv", "nv"))),
                f"Manifest geometry mismatch: {name}")
        needs_change = max(a["p"], a["q"]) > 5 or a["nu"] * a["nv"] >= 200
        actual_change = old["sha256"] != new["sha256"]
        require(actual_change == needs_change, f"Unexpected mutation or missed migration: {name}")
        require(max(b["p"], b["q"]) <= 5, f"Degree exceeded: {name}")
        require(a["weights"] == b["weights"], f"Rationality changed: {name}")
        if needs_change and int(new["metadata"][6]):
            require(b["p"] in (3, 5) and b["q"] in (3, 5) and not b["weights"],
                    f"Refit is not a cubic/quintic polynomial: {name}")
        lo, hi = ((50, 79), (80, 109), (110, 139), (140, 169), (170, 199))[(id_ % 100) // 20]
        require(lo <= b["nu"] * b["nv"] <= hi, f"Control band exceeded: {name}")
        if actual_change:
            changed.append(id_)
        if (a["nu"], a["nv"]) != (b["nu"], b["nv"]):
            grids.append({"id": id_, "old": [a["nu"], a["nv"]], "new": [b["nu"], b["nv"]]})
        degrees[f'{b["p"]},{b["q"]}'] += 1
        rows.append({"file": name, "changed": int(actual_change),
                     "old_sha256": old["sha256"], "new_sha256": new["sha256"]})
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "hashes.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    digest = lambda field: hashlib.sha256("".join(
        f'{r["file"]} {r[field]}\n' for r in rows).encode()).hexdigest()
    summary = {"policy": "dataset500-degree5-v2", "cases": 500, "max_degree": 5,
               "controls_exclusive_upper": 200, "changed": len(changed),
               "unchanged": 500 - len(changed), "changed_ids": changed,
               "grid_changes": grids, "degrees": dict(sorted(degrees.items())),
               "all_metadata_analytic_parameters_and_sources_preserved": True,
               "original_dataset_sha256": digest("old_sha256"),
               "migrated_dataset_sha256": digest("new_sha256"),
               "digest_format": "SHA256 of sorted 'case-NNN.b2a SHA256\\n' UTF-8 lines"}
    (args.output / "identity-audit.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
