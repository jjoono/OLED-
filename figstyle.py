"""Common figure style for the Optics Express submission.

Figures are drawn at the final size of the Optica single-column layout (text width 5.25 in), with sans-serif lettering of
7-9 pt, boxed axes with inward ticks, panel labels (a), (b), ... and no descriptive titles (descriptions go in the caption).
Usage:  import figstyle as FS; FS.apply(); fig = plt.figure(figsize=(FS.W, height))
"""
import matplotlib as mpl
from matplotlib.lines import Line2D

W = 5.25          # text width, in
INK = '#000000'; INK2 = '#4d4d4d'
SMALL = 7.0       # smallest lettering used inside panels


def apply():
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Liberation Sans', 'Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8, 'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
        'legend.fontsize': 7, 'legend.title_fontsize': 7, 'legend.frameon': False, 'legend.handlelength': 1.6,
        'legend.handletextpad': 0.5, 'legend.labelspacing': 0.35, 'legend.borderaxespad': 0.4,
        'axes.linewidth': 0.8, 'axes.spines.top': True, 'axes.spines.right': True, 'axes.labelpad': 3,
        'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
        'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.minor.size': 1.7, 'ytick.minor.size': 1.7,
        'xtick.major.width': 0.7, 'ytick.major.width': 0.7, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5,
        'lines.linewidth': 1.3, 'lines.markersize': 4.5,
        'mathtext.fontset': 'custom', 'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
        'mathtext.bf': 'Liberation Sans:bold', 'mathtext.fallback': 'stixsans',
        'savefig.dpi': 600, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 150,
    })


def label(ax, letter, text='', x=0.0, y=1.035, **kw):
    """Panel label '(a)' above the top-left corner of the axes, optionally followed by a short descriptor."""
    ax.text(x, y, f'({letter})' + (f' {text}' if text else ''), transform=ax.transAxes, fontsize=9, va='bottom', ha='left', **kw)


def handles(classes):
    """Legend handles for a dict {key: (color, marker, label)}."""
    return [Line2D([], [], ls='none', marker=mk, color=col, ms=5, label=lab[0].upper() + lab[1:]) for col, mk, lab in classes.values()]


def save(fig, path_no_ext):
    fig.savefig(path_no_ext + '.png', dpi=600, bbox_inches='tight', pad_inches=0.03)
    fig.savefig(path_no_ext + '.pdf', bbox_inches='tight', pad_inches=0.03)
