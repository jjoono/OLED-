"""Figure for tl_bphen10_bspline.py: n,k from the visible fit + B-spline expansion
(free spline) against the material's own Gen-Osc, for 2_23 and the two controls."""
import os as _os
ELLIPS_OUT = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_genosc as GO, tl_ce_thickness as T
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

W = _os.path.join(ELLIPS_OUT, '')
INK, INK2, MUTED, GRID, AXIS, SURF = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
C_SPL, C_REF = '#2a78d6', '#eb6834'
TITLE = {'2_23': '2_23 - nominal Bphen 10 nm', 'Src10': 'Src10 - NPB ~10 nm (control)', '2_21': '2_21 - Bphen ~47 nm (control)'}

def main():
    res = {r['piece']: r for r in json.load(open(W + 'bphen10_bspline.json'))}
    fig, ax = plt.subplots(2, 3, figsize=(11, 5.6), sharex=True, facecolor=SURF)
    for j, pc in enumerate(('2_23', 'Src10', '2_21')):
        r = res[pc]; e = r['free']; wl = np.array(e['wl'])
        Nr = np.sqrt(GO.eps_genosc(np.array(T.disp(r['mat'])['p']), wl).astype(complex))
        for i, (ys, yr, lab) in enumerate(((e['n'], Nr.real, 'n'), (e['k'], Nr.imag, 'k'))):
            a = ax[i, j]; a.set_facecolor(SURF)
            a.plot(wl, yr, color=C_REF, lw=2, label='%s Gen-Osc (tooling fits)' % r['mat'])
            a.plot(wl, ys, color=C_SPL, lw=2, label='B-spline, visible fit then expanded')
            a.grid(color=GRID, lw=0.6); a.tick_params(colors=MUTED, labelsize=8)
            for s in a.spines.values(): s.set_color(AXIS)
            a.set_ylabel(lab, color=INK2, fontsize=9)
            if i == 0:
                a.set_title('%s\nd = %.1f nm (visible fit), MSE %.1f' % (TITLE[pc], e['d'], e['mse']), color=INK, fontsize=9, loc='left')
            if pc == '2_23' and i == 1:
                for x, t, dx, ha in ((365, 'Si E1', 6, 'left'), (292, 'Si E2', -6, 'right')):
                    a.axvline(x, color=MUTED, lw=0.8, ls=':')
                    a.text(x + dx, 0.66, t, color=MUTED, fontsize=7, ha=ha)
        ax[1, j].set_xlabel('wavelength (nm)', color=INK2, fontsize=9)
    ax[0, 0].set_ylim(1.0, 2.3); ax[0, 1].set_ylim(1.0, 2.3); ax[0, 2].set_ylim(1.0, 2.3)
    for a in ax[1]: a.set_ylim(-0.02, 0.7)
    h, l = ax[0, 0].get_legend_handles_labels()
    fig.legend([h[1], h[0]], ['B-spline: visible fit (420-1000 nm), then expanded to 245 nm', 'reference: the material\'s Gen-Osc from the tooling fits'],
               loc='lower center', ncol=2, frameon=False, fontsize=8.5, labelcolor=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = W + 'bphen10_bspline_nk.png'; fig.savefig(out, dpi=200, facecolor=SURF); print(out)

if __name__ == '__main__':
    main()
