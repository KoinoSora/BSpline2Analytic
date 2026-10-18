# 文件职责：测试代码、案例与材料的清单生成器。
# 读取现有测试证据形成可追溯清单，不以文档统计替代算法回归。

"""只读生成测试输入和源码清单，不运行几何算法。"""
from pathlib import Path
import collections
import csv
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    paths = sorted(p for p in tracked if p.endswith(".b2a"))
    manifest_by_family = {
        "iso500": "tests/data/dataset500/manifest.csv",
        "business200": "tests/data/business-extreme200/manifest.csv",
        "historical_extreme200": "tests/data/extreme200/manifest.csv",
        "retired_general600": "tests/data/dataset500/manifest-general600.csv",
        "retired_extra60": "tests/data/dataset500/manifest-general-extra60.csv",
        "retired_extreme320": "tests/data/general-extreme320/manifest.csv",
    }
    maps = {k: {int(r["id"]): r for r in read_csv(ROOT / v)} for k, v in manifest_by_family.items()}
    v2 = json.loads((ROOT / "tests/data/general-rational-v2/dataset.json").read_text())
    vmap = {(r["suite"], int(r["id"])): r for r in v2["cases"]}
    rows = []
    for rel in paths:
        path = ROOT / rel
        if "/general-rational-v2/" in rel:
            suite = path.parent.name
            case_id = int(path.stem.split("-")[-1])
            family, state = "v2_" + suite, "current"
            metadata = vmap[suite, case_id]
            assert sha(path) == metadata["sha256"], rel
            generator = "tests/GeneralRationalDataset.cpp;tests/GeneralRationalReference.cpp"
        elif rel.startswith("tests/data/"):
            case_id = int(path.stem.split("-")[-1])
            if path.parent.name == "dataset500":
                family = "retired_extra60" if path.name.startswith("general-extra-") else "retired_general600" if path.name.startswith("general-") else "iso500"
            else:
                family = {"business-extreme200": "business200", "extreme200": "historical_extreme200", "general-extreme320": "retired_extreme320"}[path.parent.name]
            metadata = maps[family][case_id]
            state = "archived_no_algorithm_runs" if family.startswith("retired_") else "current_specialized_only" if family != "iso500" else "current"
            generator = "tests/" + {"iso500": "Dataset500.cpp", "business200": "BusinessExtremeDataset.cpp", "historical_extreme200": "ExtremeDataset.cpp", "retired_general600": "GeneralDataset.cpp", "retired_extra60": "GeneralDataset.cpp", "retired_extreme320": "GeneralExtremeDataset.cpp"}[family]
        else:
            assert rel == "docs/testing/introduction-cases-2026-10-16/torus-triangle.b2a", rel
            family, case_id, state = "introduction_triangle", "torus-triangle", "current_demonstration"
            metadata, generator = {}, "tests/IntroductionCases.cpp"
        rows.append(dict(family=family, id=case_id, status=state, path=rel, sha256=sha(path),
                         generator=generator, metadata_json=json.dumps(metadata, ensure_ascii=False, sort_keys=True)))
    assert len(rows) == 2866, f"Fixture population changed: {len(rows)}; review documentation quotas."
    write_csv("case-inventory.csv", rows)
    sources = [dict(path=p, sha256=sha(ROOT / p)) for p in tracked
               if p.startswith("tests/") and Path(p).suffix in (".cpp", ".h", ".py")]
    write_csv("test-source-inventory.csv", sorted(sources, key=lambda r: r["path"]))
    summary = dict(scope="tracked .b2a fixtures; identities, not independent geometric shapes",
                   total=len(rows), families=dict(sorted(collections.Counter(r["family"] for r in rows).items())),
                   test_source_files=len(sources),
                   inventory_sha256=sha(OUT / "case-inventory.csv"),
                   source_inventory_sha256=sha(OUT / "test-source-inventory.csv"),
                   geometry_executed=False)
    (OUT / "inventory-summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
