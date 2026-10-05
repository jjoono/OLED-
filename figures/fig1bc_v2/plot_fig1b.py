"""Redesigned Fig. 1b: one y-axis (fraction of substrate light), per-pass bars + cumulative lines, direct labels."""
import os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
P = 0.4
CASES = [(0.02, 'low loss', '#2a78d6'), (0.045, 'Ag', '#1baf7a'), (0.15, 'Al', '#eb6834')]
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.linewidth': 0.6,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5})
K = np.arange(0, 11); N = 30
fig, ax = plt.subplots(figsize=(3.5, 2.4))   # 89 mm single column
bw = 0.24
for i, (a, lab, col) in enumerate(CASES):
    q = (1 - P) * (1 - a)
    per = P * q ** K
    cum = np.cumsum(P * q ** np.arange(N))
    ax.bar(K + (i - 1) * (bw + 0.02), per, bw, color=col, alpha=0.45, linewidth=0, zorder=2)
    ax.plot(K, cum[:len(K)], color=col, lw=1.6, zorder=3, solid_capstyle='round')
    lim = P / (P + (1 - P) * a)
    ylab = (1.00, 0.89, 0.78)[i]
    ax.annotate(f"{a:<5}  {lab}  → {lim:.2f}", xy=(10, cum[10]), xytext=(10.7, ylab), textcoords='data',
                fontsize=6.3, color='#0b0b0b', va='center', ha='left', annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color=col, lw=1.0, shrinkA=0, shrinkB=1))
ax.set_xlim(-0.6, 10.3); ax.set_ylim(0, 1.0)
ax.set_xticks(range(0, 11, 2)); ax.set_yticks(np.arange(0, 1.01, 0.2))
ax.set_xlabel('Number of round trips'); ax.set_ylabel('Fraction of substrate light')
ax.grid(axis='y', color='#e6e5e0', lw=0.5, zorder=0); ax.set_axisbelow(True)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.text(10.75, 1.065, 'A′', fontsize=6.3, color='#52514e', ha='left', clip_on=False)
ax.text(4.3, 0.56, f'$p$ = {P}', fontsize=7, color='#0b0b0b')
ax.text(4.3, 0.47, 'Lines: cumulative extraction', fontsize=6.3, color='#52514e')
ax.text(4.3, 0.40, 'Bars: extracted per round trip', fontsize=6.3, color='#52514e')
fig.subplots_adjust(left=0.13, right=0.70, bottom=0.17, top=0.93)
for ext in ('png', 'pdf', 'svg'): fig.savefig(os.path.join(HERE, f'fig1b_redesign.{ext}'), dpi=400, bbox_inches='tight', pad_inches=0.04)
