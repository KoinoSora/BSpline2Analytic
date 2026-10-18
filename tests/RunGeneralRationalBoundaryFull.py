# 文件职责：当前 V2 一般兜底的完整串行回归调度。
# 逐例记录执行失败、超时和审计结果，保存可续跑证据；不得调用旧验证器作对照。

"""完整覆盖当前有理边界验证器，不调用任何废弃验证器。
逐例串行进程避免并发基准争用；每例三次全新上下文运行，来源证书、输出审计与已知母面诊断均在算法计时外。
中断及失败明确记录，不能填零冒充测量。"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

SIZES = dict(main=600, extra=60, extreme=320, synthetic=5)


def read(path):
    if not path.exists():
        return []
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('exe', 'sdk', 'data', 'output'):
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--fixed-chart', action='store_true', help='Separate supplied-mother/chart diagnostics; never classified timing')
    p.add_argument('--runs', type=int, default=3)
    p.add_argument('--timeout', type=float, default=60)
    p.add_argument('--max-case-ms', type=float, default=50,
                   help='Shared classification + validation wall budget; 0 disables only for explicit diagnostics')
    p.add_argument('--cases', nargs='*', help='Optional smoke-test keys, e.g. main:100 extreme:0')
    p.add_argument('--resume', action='store_true')
    args = p.parse_args()
    if args.runs < 1 or not 0 < args.timeout <= 300:
        p.error('Invalid runs/timeout')
    if not math.isfinite(args.max_case_ms) or not 0 <= args.max_case_ms <= 60000:
        p.error('Invalid max-case-ms')
    exe, sdk, data, output = (getattr(args, name).resolve() for name in ('exe', 'sdk', 'data', 'output'))
    from RationalDatasetSupport import require_v2_data
    require_v2_data(data)
    cases = [(s, i) for s, size in SIZES.items() for i in range(size)]
    if args.cases:
        selected = [tuple(key.split(':')) for key in args.cases]
        cases = [(s, int(i)) for s, i in selected]
        if len(set(cases)) != len(cases) or any(s not in SIZES or not 0 <= i < SIZES[s] for s, i in cases):
            p.error('Invalid/duplicate case key')
    metadata = {}
    for suite in SIZES:
        directory, manifest = data / suite, 'manifest.csv'
        entries = read(directory / manifest)
        if [int(r['id']) for r in entries] != list(range(SIZES[suite])):
            raise ValueError(f'Incomplete/out-of-order manifest: {suite}')
        for row in entries:
            metadata[(suite, int(row['id']))] = row
    inputs = {}
    for suite, i in cases:
        inputs[f'{suite}:{i}'] = sha(data / suite / metadata[(suite, i)]['file'])
    config = dict(validator='rational_boundary', dataset_contract='rational-domain-v2', fixed_chart=args.fixed_chart, exe_sha256=sha(exe), inputs_sha256=inputs,
                  expectations={f"{s}:{i}": metadata[(s, i)]["expectation"] for s, i in cases},
                  runs=args.runs, timeout_seconds=args.timeout, max_case_ms=args.max_case_ms,
                  cancellation_reserve_ms=min(5.,args.max_case_ms*.1), cases=cases,
                  cpu_affinity_fixed=False, exclusive_machine=False, parallel_processes=1)
    output.mkdir(parents=True, exist_ok=True)
    stamp = output / 'run.json'
    if args.resume:
        if not stamp.exists() or json.loads(stamp.read_text()) != json.loads(json.dumps(config)):
            raise ValueError('Resume requires identical executable, inputs and run configuration')
    elif stamp.exists():
        raise ValueError('Output already used; choose a new output or explicit --resume')
    stamp.write_text(json.dumps(config, indent=2), encoding='utf-8')
    parts = output / 'cases'
    parts.mkdir(exist_ok=True)
    env = dict(os.environ)
    key = next(k for k in env if k.lower() == 'path')
    env[key] = str(sdk / 'x64-win/bin') + os.pathsep + env[key]
    all_rows, all_attempts, known, executions = [], [], [], []
    for index, (suite, case) in enumerate(cases):
        target = parts / f'{suite}-{case:03}.csv'
        log, done = target.with_suffix('.log'), target.with_suffix('.done.json')
        if args.resume and done.exists():
            execution = json.loads(done.read_text())
        else:
            target.unlink(missing_ok=True)
            start = time.perf_counter()
            with log.open('wb') as capture:
                try:
                    result = subprocess.run([str(exe), '--rational-fixed-case' if args.fixed_chart else '--rational-full-case', str(data),
                                             str(target), suite, str(case), str(args.runs), '--rational-v2-only', str(args.max_case_ms)],
                                            cwd=sdk, env=env, stdout=capture, stderr=subprocess.STDOUT,
                                            timeout=args.timeout, creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
                    code, timeout = result.returncode, False
                except subprocess.TimeoutExpired:
                    code, timeout = None, True
            content = log.read_text(encoding='utf-8', errors='replace')
            markers = [line for line in content.splitlines() if line.startswith(('BEGIN_', 'END_'))]
            execution = dict(suite=suite, id=case, exit_code=code, timeout=timeout,
                             wall_seconds=time.perf_counter()-start, last_marker=markers[-1] if markers else '',
                             log_tail=content[-1500:] if timeout or code else '')
            done.write_text(json.dumps(execution, indent=2), encoding='utf-8')
        executions.append(execution)
        rows = read(target)
        if any(r.get('suite') != suite or r.get('id') != str(case) or r.get('mode') != 'rational' for r in rows):
            raise ValueError(f'Invalid output identity/mode: {target}')
        if [int(r['run']) for r in rows] != list(range(len(rows))) or len(rows) > args.runs:
            raise ValueError(f'Invalid output run sequence: {target}')
        for row in rows:
            if float(row.get('max_case_ms', 'nan')) != args.max_case_ms:
                raise ValueError('Executable did not apply the requested case budget')
            row['input_sha256'] = inputs.get(f'{suite}:{case}', 'constructed')
            row['execution_status'] = 'timeout' if execution['timeout'] else 'process_error' if execution['exit_code'] else 'completed'
            row['watchdog_seconds'] = execution['wall_seconds']
        for run in range(len(rows), args.runs):
            meta = metadata.get((suite, case), {})
            rows.append(dict(suite=suite, id=str(case), run=str(run), mode='rational',
                             group=meta.get('group', 'nurbs_chart_pressure'), source_type=meta.get('type', ''),
                             rational=meta.get('rational', ''), degree_u=meta.get('degree_u', ''), degree_v=meta.get('degree_v', ''),
                             controls=meta.get('controls', ''), status='not_measured', accepted='', audit_status='not_available',
                             total_ms='nan', message=execution['last_marker'],
                             execution_status='timeout' if execution['timeout'] else 'process_error',
                             input_sha256=inputs.get(f'{suite}:{case}', 'constructed'), watchdog_seconds=execution['wall_seconds']))
        all_rows.extend(rows)
        all_attempts.extend(read(target.with_name(target.stem + '-attempts.csv')))
        known.extend(read(target.with_name(target.stem + '-known.csv')))
        # 逐例写入检查点，包括执行失败及部分计时。
        write(output / 'results.csv', all_rows)
        if index % 10 == 0 or index + 1 == len(cases):
            passed = sum(r.get('accepted') == '1' for r in rows)
            print(f'{index+1}/{len(cases)} {suite}:{case}: accepted={passed}/{args.runs}, status={rows[-1]["status"]}, wall={execution["wall_seconds"]:.2f}s', flush=True)
    write(output / 'attempts.csv', all_attempts)
    write(output / 'known-mother.csv', known)
    (output / 'executions.json').write_text(json.dumps(executions, indent=2), encoding='utf-8')
    errors = [x for x in executions if x['timeout'] or x['exit_code']]
    print(f'Completed {len(cases)} cases; {len(all_rows)} rows; {len(errors)} execution errors', flush=True)
    overruns = [r for r in all_rows if r.get('within_time_limit') == '0']
    print(f'Wall-budget violations: {len(overruns)} runs; cooperative cancellations: '
          f'{sum(r.get("time_budget_exceeded") == "1" for r in all_rows)} runs', flush=True)
    return int(bool(errors or overruns))


if __name__ == '__main__':
    raise SystemExit(main())
