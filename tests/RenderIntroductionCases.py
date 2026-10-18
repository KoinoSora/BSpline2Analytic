# 文件职责：将实测转换导出数据渲染为文档插图。
# 源/目标保持相同相机和比例，不合成或替换 CAD 几何；文件身份用于关联真实案例。

"""渲染真实串联转换证据，不合成替代 CAD 几何。
输入网格及输出边界采样来自 IntroductionCases.cpp；目标网格由恢复的母面、参数图和多边形求值，不是 BRep 面对象。
热力图表示采样母面距离，不是连续域误差证书。"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import font_manager
import numpy as np


def rows(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--evidence', type=Path, required=True)
    ap.add_argument('--figures', type=Path, required=True)
    args = ap.parse_args()
    repo = Path(__file__).resolve().parents[1]
    args.figures.mkdir(parents=True, exist_ok=True)
    font = Path('C:/Windows/Fonts/msyh.ttc')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family'] = font_manager.FontProperties(fname=str(font)).get_name()
    plt.rcParams.update({'axes.unicode_minus': False, 'font.size': 10, 'svg.fonttype': 'path'})
    run_rows = rows(args.evidence / 'runs.csv')
    names = {'sphere-317': '球面极点 · case-317', 'cone-217': '圆锥顶点 · case-217',
             'torus-triangle': '规则圆环面 · 三角形裁剪与坍缩边'}
    summaries = []
    for name, title in names.items():
        runs = [r for r in run_rows if r['case'] == name]
        assert len(runs) == 3
        data = rows(args.evidence / (name + '-geometry.csv'))
        xyz = lambda part: np.array([[float(r[k]) for k in ('x', 'y', 'z')] for r in part])
        source = [r for r in data if r['kind'] == 'source']
        n = max(int(r['i']) for r in source) + 1
        source_grid = xyz(source).reshape(n, n, 3)
        target_grid = xyz([r for r in data if r['kind'] == 'target']).reshape(n, n, 3)
        assert np.isfinite(source_grid).all() and np.isfinite(target_grid).all()
        errors = np.array([float(r['mother_distance']) for r in source]).reshape(n, n)
        summary = dict(case=name, runs=3, all_accepted=all(r['accepted'] == '1' for r in runs),
                       all_audited=all(r['audit'] == '1' for r in runs),
                       min_ms=min(float(r['total_ms']) for r in runs),
                       mean_ms=float(np.mean([float(r['total_ms']) for r in runs])),
                       max_ms=max(float(r['total_ms']) for r in runs),
                       boundary_bound=max(float(r['boundary_bound']) for r in runs),
                       tangent_bound_deg=max(float(r['tangent_bound_deg']) for r in runs),
                       plot_sample_max_mm=float(errors.max()), degree_u=int(runs[0]['degree_u']),
                       degree_v=int(runs[0]['degree_v']), controls=int(runs[0]['controls']),
                       collapsed_edges=int(runs[0]['collapsed_edges']),
                       micro_exemptions=int(runs[0]['micro_exemptions']))
        summary['under_10ms'] = summary['max_ms'] < 10
        summary['within_50ms'] = summary['max_ms'] <= 50
        summaries.append(summary)
        fig = plt.figure(figsize=(12.6, 4.8), facecolor='white')
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, .82], left=.025, right=.96, bottom=.20, top=.80, wspace=.10)
        points = np.concatenate([source_grid.reshape(-1, 3), target_grid.reshape(-1, 3)])
        lo, hi = points.min(axis=0), points.max(axis=0)
        mid = (lo + hi) / 2
        radius = max(hi - lo) * .57
        normal = np.cross(target_grid[n//2+1,n//2]-target_grid[n//2-1,n//2],
                          target_grid[n//2,n//2+1]-target_grid[n//2,n//2-1])
        normal /= np.linalg.norm(normal)
        if normal[2] < 0:
            normal = -normal
        elev = min(75, float(np.degrees(np.arcsin(normal[2]))) + 10)
        azim = float(np.degrees(np.arctan2(normal[1], normal[0]))) + 12
        # 源面和目标面使用同一正交相机及相同几何比例。
        for k, grid in enumerate((source_grid, target_grid)):
            ax = fig.add_subplot(gs[0, k], projection='3d', proj_type='ortho', computed_zorder=False)
            ax.plot_surface(grid[:, :, 0], grid[:, :, 1], grid[:, :, 2],
                            color=('#79b9d4' if k == 0 else '#82bfaa'), alpha=1,
                            rcount=n, ccount=n, linewidth=0, antialiased=False, shade=True, zorder=1)
            if k == 0:
                ax.plot_wireframe(*np.moveaxis(grid, 2, 0), rstride=5, cstride=5,
                                  color='#386988', linewidth=.45, alpha=.65, zorder=2)
            kind = 'source_edge' if k == 0 else 'output_edge'
            for edge in range(4):
                part = xyz([r for r in data if r['kind'] == kind and int(r['edge']) == edge])
                ax.plot(*part.T, color=('#174c71' if k == 0 else '#bb5424'), linewidth=1.8, zorder=3)
            for r in data:
                if r['kind'] == 'collapsed':
                    point = [float(r[key]) for key in ('x', 'y', 'z')]
                    ax.scatter(*point, color='#c42f42', s=55, depthshade=False, zorder=20)
            ax.set(xlim=(mid[0]-radius, mid[0]+radius), ylim=(mid[1]-radius, mid[1]+radius),
                   zlim=(mid[2]-radius, mid[2]+radius))
            ax.set_box_aspect((1, 1, 1), zoom=(1.1 if name == 'torus-triangle' else 1.3)); ax.view_init(elev=elev, azim=azim)
            ax.set_axis_off()
            ax.set_title('输入样条面与参数网格' if k == 0 else '恢复母面与实际输出边界', fontsize=11, pad=0)
        ax = fig.add_subplot(gs[0, 2])
        im = ax.imshow(errors.T, origin='lower', extent=(0, 1, 0, 1), aspect='equal',
                       cmap='viridis', vmin=0, vmax=max(float(errors.max()), 1e-15), interpolation='nearest')
        ax.set_title('采样点到恢复母面的距离', fontsize=11, pad=12)
        ax.set_xlabel('归一化源参数 u'); ax.set_ylabel('归一化源参数 v')
        ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
        cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.045)
        cb.formatter.set_powerlimits((0, 0)); cb.update_ticks(); cb.set_label('mm')
        fig.suptitle(title, x=.035, y=.975, ha='left', fontsize=17, fontweight='bold', color='#15394e')
        fig.text(.035, .876, f"{summary['degree_u']} × {summary['degree_v']} 次  |  {summary['controls']} 控制顶点  |  等参 → 兜底  |  3/3 成功且审计通过  |  整体最大耗时 {summary['max_ms']:.4f} ms", color='#405560', fontsize=10)
        fig.legend(handles=[Line2D([], [], color='#174c71', label='源边界'),
                            Line2D([], [], color='#bb5424', label='实际输出边界'),
                            Line2D([], [], marker='o', color='w', markerfacecolor='#c42f42', label='整条边坍缩至此点', markersize=7)],
                   loc='lower left', bbox_to_anchor=(.025, .075), ncol=3, frameon=False, fontsize=9)
        fig.text(.035, .045, f"几何按真实比例显示；热图为 {n}×{n} 离散采样，非连续域上界。绘图与输出均在算法计时之外。", fontsize=9, color='#586570')
        fig.savefig(args.figures / (name + '.png'), dpi=180)
        fig.savefig(args.figures / (name + '.svg'))
        svg = args.figures / (name + '.svg')
        svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines()) + '\n', encoding='utf-8')
        plt.close(fig)
    (args.evidence / 'summary.json').write_text(json.dumps(summaries, ensure_ascii=False, indent=2), encoding='utf-8')
    inputs = [repo/'tests/data/dataset500/case-317.b2a', repo/'tests/data/dataset500/case-217.b2a',
              repo/'tests/data/general-rational-v2/main/main-400.b2a', args.evidence/'torus-triangle.b2a']
    source_files = sorted([*repo.glob('src/*.cpp'), *repo.glob('src/*.h'), *repo.glob('include/*.h')])
    source_files += [repo/'tests/IntroductionCases.cpp', repo/'tests/SequentialOutputAudit.h', Path(__file__).resolve(), repo/'CMakeLists.txt']
    source_files += sorted([*repo.glob('tests/*.cpp'), *repo.glob('tests/*.h')])
    exe = repo/'build-sequential-20261016/Release/BSpline2AnalyticIntroductionCases.exe'
    provenance = dict(base_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
                      platform=platform.platform(), processor=os.environ.get('PROCESSOR_IDENTIFIER'),
                      executable_sha256=sha(exe), source_sha256={str(p.relative_to(repo)):sha(p) for p in source_files},
                      input_sha256={str(p.resolve().relative_to(repo)):sha(p) for p in inputs},
                      artifact_sha256={str(p.resolve().relative_to(repo)):sha(p) for p in sorted(args.evidence.glob('*.csv'))},
                      timing='External steady_clock across complete sequential API; 3 fresh contexts per case; no per-case warmup; audits/geometry exports between runs, outside timers; one serial process; includes first call.',
                      plotting='Source: actual input evaluations; target: recovered analytic mother with recovered rational chart/corners; boundary: actual SDK output curve evaluations. Same local output frame, same camera/scale. No geometric displacement amplification.',
                      limitations='Selected-case evidence only. Not a new full regression; no adjacent-face G1 or continuous interior proof. All selected cases used zero optimization iterations.',
                      matplotlib=matplotlib.__version__, numpy=np.__version__)
    (args.evidence/'provenance.json').write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(summaries, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
