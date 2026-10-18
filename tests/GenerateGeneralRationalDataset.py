# 文件职责：V2 输入生成、来源认证与可续跑调度。
# 只生成已认证数据，不执行生产转换或验证测试；续跑校验现有文件身份。

"""仅生成并独立认证 V2 输入，不运行转换或生产验证测试。"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
from RationalDatasetSupport import REPRESENTATION_POLICY, require_representation

SIZES = dict(main=600, extra=60, extreme=320, synthetic=5)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('exe', 'sdk', 'legacy-main', 'legacy-extreme', 'output'):
        p.add_argument('--' + key, type=Path, required=True)
    p.add_argument('--suites', nargs='+', choices=SIZES, default=list(SIZES))
    p.add_argument('--resume', action='store_true')
    p.add_argument('--audit', action='store_true', help='Re-certify stored geometry and record SDK sample/size diagnostics only')
    a = p.parse_args()
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    logs = root / 'generation-logs'
    logs.mkdir(exist_ok=True)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    env[key] = str(a.sdk.resolve() / 'x64-win/bin') + os.pathsep + env[key]
    generated = False
    for suite in a.suites:
        existing = [i for i in range(SIZES[suite]) if (root / suite / f'{suite}-{i:03}.b2a').exists()]
        if existing and not a.resume:
            p.error('Existing V2 inputs: explicit --resume required')
        missing = sorted(set(range(SIZES[suite])) - set(existing))
        if not missing:
            continue
        # 批次失败时只发布已独立认证的文件；续跑从首个缺失文件开始，
        # 不能静默覆盖已有输入数据。
        first = missing[0]
        while first < SIZES[suite]:
            if first in existing:
                first += 1
                continue
            end = first + 1
            while end < SIZES[suite] and end not in existing:
                end += 1
            log = logs / f'{suite}-{first:03}-{end:03}.log'
            cmd = [str(a.exe.resolve()), '--generate', str(a.legacy_main.resolve()),
                   str(a.legacy_extreme.resolve()), str(root), suite, str(first), str(end)]
            print(f'GENERATING {suite} [{first}, {end})', flush=True)
            with log.open('wb') as out:
                result = subprocess.run(cmd, cwd=a.sdk.resolve(), env=env, stdout=out,
                                        stderr=subprocess.STDOUT,
                                        creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            print(log.read_text(encoding='utf-8', errors='replace')[-2500:], flush=True)
            if result.returncode:
                raise RuntimeError(f'Data generation failed; see {log}')
            generated = True
            first = end
    if a.audit:
        for suite in a.suites:
            log = logs / f'audit-{suite}.log'
            cmd = [str(a.exe.resolve()), '--audit', str(root), suite, '0', str(SIZES[suite])]
            print('SOURCE AUDIT ' + suite, flush=True)
            with log.open('wb') as out:
                result = subprocess.run(cmd, cwd=a.sdk.resolve(), env=env, stdout=out,
                                        stderr=subprocess.STDOUT,
                                        creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            if result.returncode:
                raise RuntimeError(f'Source audit failed; see {log}')
    records = []
    for suite, size in SIZES.items():
        manifest = root / suite / 'manifest.csv'
        if not manifest.exists():
            print('Partial generation; dataset.json is not published', flush=True)
            return
        with manifest.open(encoding='utf-8', newline='') as stream:
            rows = list(csv.DictReader(stream))
        if [int(r['id']) for r in rows] != list(range(size)):
            raise ValueError('Incomplete manifest: ' + suite)
        for row in rows:
            path = root / suite / row['file']
            if path.read_text().splitlines()[0] != 'B2A_RATIONAL_DOMAIN_V2':
                raise ValueError('Retired input in V2 corpus')
            actual = require_representation(path.read_bytes())
            if actual != (int(row['degree_u']), int(row['degree_v']), int(row['controls'])):
                raise ValueError('Manifest representation mismatch')
            records.append(dict(suite=suite, **row, sha256=sha(path)))
    previous = json.loads((root/'dataset.json').read_text()) if (root/'dataset.json').exists() else {}
    stamp = dict(contract='rational-domain-v2', cases=records, sizes=SIZES,
                 representation_policy=REPRESENTATION_POLICY,
                 last_generation_exe_sha256=sha(a.exe) if generated else previous.get('last_generation_exe_sha256', previous.get('generation_exe_sha256')),
                 source_certification='independent_rational_domain',
                 fallback_tests_run=False)
    if a.audit and set(a.suites) == set(SIZES):
        stamp['qualification_exe_sha256'] = sha(a.exe)
        stamp['qualification_csv_sha256'] = {suite: sha(root/suite/f'audit-0-{size}.csv') for suite, size in SIZES.items()}
    (root / 'dataset.json').write_text(json.dumps(stamp, indent=2, ensure_ascii=False), encoding='utf-8')
    print(f'Published {len(records)} certified V2 inputs; no fallback tests run', flush=True)


if __name__ == '__main__':
    main()
