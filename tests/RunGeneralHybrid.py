# 文件职责：已废弃混合验证器的历史逐例调度器。
# 主函数主动退出，不得调用 archived_main 绕过禁测。

"""旧采样/连续混合验证器的串行逐例看门狗，当前禁用。
历史实现只用真实分类候选，真值不进入生产代码；来源/输出证书、加载与进程启动均在分类加验证计时外，不调用优化器。
不得通过 archived_main 绕过禁测。"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def summarize(rows, failures):
    groups = {}
    for row in rows:
        groups.setdefault(row['group'], []).append(row)
    result = {}
    for group, items in groups.items():
        per_case = {}
        for row in items:
            per_case.setdefault((row['extra'], row['id']), []).append(row)
        values = sorted(statistics.median(float(r['pipeline_ms']) for r in rs)
                        for rs in per_case.values())
        def p(q):
            return values[max(0, math.ceil(q * len(values)) - 1)]
        result[group] = {
            'completed_cases': len(per_case), 'runs': len(items),
            'accepted_cases_all_runs': sum(all(r['accepted'] == '1' for r in rs) for rs in per_case.values()),
            'audit_failures': sum(r['accepted'] == '1' and r['output_audit'] != '1' for r in items),
            'pipeline_case_median_mean_ms': statistics.mean(values),
            'pipeline_case_median_p50_ms': p(.50), 'pipeline_case_median_p95_ms': p(.95),
            'pipeline_case_median_p99_ms': p(.99), 'pipeline_case_median_max_ms': max(values),
            'cases_median_le_10ms': sum(v <= 10 for v in values),
            'accepted_cases_all_runs_le_10ms': sum(all(r['accepted'] == '1' and float(r['pipeline_ms']) <= 10 for r in rs)
                                                 for rs in per_case.values()),
            'unstable_outcome_cases': sum(len({(r['accepted'], r['status'], r['output_type'], r['output_audit']) for r in rs}) > 1
                                         for rs in per_case.values()),
            'mean_stages_ms': {key: statistics.mean(float(r[key]) for r in items)
                               for key in ('classify_ms', 'sampling_ms', 'recovery_ms', 'prepare_ms',
                                           'boundary_ms', 'validate_ms')},
        }
    return {'groups': result, 'watchdog_or_process_failures': failures}


def archived_main():
    p = argparse.ArgumentParser()
    p.add_argument('--exe', required=True, type=Path)
    p.add_argument('--sdk', required=True, type=Path)
    p.add_argument('--data', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--first', type=int, default=0)
    p.add_argument('--last', type=int, default=660)
    p.add_argument('--runs', type=int, default=3)
    p.add_argument('--timeout', type=float, default=30)
    a = p.parse_args()
    a.exe, a.sdk, a.data, a.output = [x.resolve() for x in (a.exe, a.sdk, a.data, a.output)]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    work = a.output.parent / (a.output.stem + '-cases')
    work.mkdir(exist_ok=True)
    env = dict(os.environ)
    path_key = next(k for k in env if k.lower() == 'path')
    env[path_key] = str(a.sdk / 'x64-win/bin') + os.pathsep + env[path_key]
    rows, failures, input_hashes = [], [], {}
    started = time.monotonic()
    for index in range(a.first, a.last):
        extra, case_id = int(index >= 600), index - 600 if index >= 600 else index
        name = f'general-extra-{case_id:03}.b2a' if extra else f'general-{case_id:03}.b2a'
        input_hashes[name] = digest(a.data / name)
        target = work / f'{index:03}.csv'
        if target.exists():
            target.unlink()  # 只处理当前生成的逐例输出，不能复用过期记录。
        command = [str(a.exe), '--case', str(a.data), str(target), str(case_id), str(extra), str(a.runs)]
        try:
            run = subprocess.run(command, cwd=a.sdk, env=env, capture_output=True, text=True,
                                 encoding='utf-8', errors='replace', timeout=a.timeout)
            (work / f'{index:03}.log').write_text(run.stdout + run.stderr, encoding='utf-8')
            if run.returncode:
                failures.append({'index': index, 'reason': 'process_failure', 'code': run.returncode})
        except subprocess.TimeoutExpired as error:
            partial = error.stdout or b''
            if isinstance(partial, bytes):
                partial = partial.decode('utf-8', errors='replace')
            (work / f'{index:03}.log').write_text(partial, encoding='utf-8')
            failures.append({'index': index, 'reason': 'timeout', 'seconds': a.timeout,
                             'last_stage': partial.splitlines()[-1] if partial else 'startup'})
        current = list(csv.DictReader(target.open(encoding='utf-8'))) if target.exists() else []
        # 保留部分测量证据，但不能计为完整案例。
        if len(current) == a.runs:
            rows.extend(current)
        elif not any(f['index'] == index for f in failures):
            failures.append({'index': index, 'reason': 'incomplete_rows', 'rows': len(current)})
        if rows:
            with a.output.open('w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)
        if (index-a.first) % 10 == 0 or index+1 == a.last:
            print(f'Completed {index+1-a.first}/{a.last-a.first}; rows={len(rows)} failures={len(failures)}', flush=True)
    summary = summarize(rows, failures)
    summary.update(scope='Classified candidates + standalone hybrid validation; no optimizer or legacy validator',
                   repeats=a.runs, requested_cases=a.last-a.first, timeout_seconds=a.timeout,
                   wall_seconds=time.monotonic()-started, exe_sha256=digest(a.exe),
                   source_data=str(a.data), configuration='Release /fp:strict',
                   cpu_affinity_fixed=False, exclusive_machine=False)
    repo = Path(__file__).resolve().parent.parent
    sources = ['include/GeneralHybridValidation.h', 'src/GeneralHybridInternal.h',
               'src/GeneralHybridValidation.cpp', 'tests/GeneralHybridTests.cpp', 'tests/RunGeneralHybrid.py']
    a.output.with_suffix('.provenance.json').write_text(json.dumps({
        'inputs_sha256': input_hashes,
        'sources_sha256': {name: digest(repo / name) for name in sources},
        'exe_sha256': summary['exe_sha256'],
    }, indent=2), encoding='utf-8')
    a.output.with_suffix('.summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    return int(bool(failures) or any(g['audit_failures'] for g in summary['groups'].values()))


def main():
    raise SystemExit('ValidateGeneralHybrid is retired. Historical runner retained; execution disabled. Use RunGeneralRationalBoundarySelected.py.')


if __name__ == '__main__':
    raise SystemExit(main())
