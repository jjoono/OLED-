"""ZnS on Si: n,k table and figure from zns_fits.json -> $ELLIPS_OUT/zns/
  ZnS_nk.csv        wavelength_nm, n, k (depth-averaged, joint graded fit), n at the top and at the bottom of the film
  zns_summary.png   Psi / Delta of both pieces with the model, the thickness scan, n and k"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
_H = _os.path.dirname(_os.path.abspath(__file__)); sys.path.insert(0, _H); sys.path.insert(0, _os.path.join(_H, '../grahil_perov_261124'))
import zn_load as L, zn_fit as Z, gp_model as P
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

W = _os.path.join(ELLIPS_OUT, ''); OUT = _os.path.join(ELLIPS_OUT, 'zns', '')
INK, INK2, MUTED, GRID, AXIS, SURF = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'


def style(a):
    a.set_facecolor(SURF); a.grid(color=GRID, lw=0.6); a.tick_params(colors=MUTED, labelsize=8)
    for s in a.spines.values(): s.set_color(AXIS)


def main():
    _os.makedirs(OUT, exist_ok=True); J = json.load(open(W + 'zns_fits.json'))
    pick = lambda o: [r for r in J['fits'] if len(r['sheets']) == 2 and tuple(r['opts']) == o][0]
    base, gr = pick(()), pick(('graded',)); v = gr['struct']; g = v['g']
    wl = np.arange(205.0, 1000.5, 1.0); N = Z.nk_of(gr, wl); chi = N**2 - 1
    Nt, Nb = np.sqrt(1 + (1 - g) * chi), np.sqrt(1 + (1 + g) * chi)
    np.savetxt(OUT + 'ZnS_nk.csv', np.c_[wl, N.real, N.imag, Nt.real, Nb.real], delimiter=',', fmt='%.5f', comments='', header='wavelength_nm,n,k,n_top,n_bottom')
    lab = dict(color=INK2, fontsize=8.5); ttl = dict(color=INK, fontsize=9.5, loc='left')
    fig = plt.figure(figsize=(13, 9), facecolor=SURF); gs = fig.add_gridspec(3, 4, hspace=0.42, wspace=0.36, left=0.06, right=0.985, top=0.93, bottom=0.07)
    S = Z.setup(list(L.SHEETS), Z.LO, Z.HI, 1); Nw = Z.nk_of(gr, S['wl'])
    for i, sh in enumerate(L.SHEETS):
        m = Z.stack(S, Nw, v['d%d' % (i + 1)], v['d_r%d' % (i + 1)], g); ps, dl = P.psidelta(*m); d = S['D'][i]
        for j, (key, mod, yl) in enumerate((('psi', ps, 'Ψ (deg)'), ('del', dl, 'Δ (deg)'))):
            a = fig.add_subplot(gs[j, 2 * i:2 * i + 2]); style(a)
            for k, c in enumerate((C1, C2, C3)):
                y = d[key][:, k] % 360 if key == 'del' else d[key][:, k]; mm = mod[:, k].copy(); y = y.copy()
                for arr in (y, mm): arr[np.r_[False, np.abs(np.diff(arr)) > 180]] = np.nan
                a.plot(S['wl'], y, color=c, lw=2.8, alpha=0.35); a.plot(S['wl'], mm, color=c, lw=1.1, label='%g°' % L.ANG[k])
            a.set_ylabel(yl, **lab)
            if j == 0:
                a.set_title('%s: ZnS %.1f nm + roughness %.1f nm, MSE %.1f   (thick line measured, thin line model)' % (sh, v['d%d' % (i + 1)], v['d_r%d' % (i + 1)], gr['mse'][i]), **ttl); a.legend(frameon=False, fontsize=7.5, labelcolor=INK2, ncol=3, loc='upper right')
            else: a.set_xlabel('wavelength (nm)', **lab)
    a = fig.add_subplot(gs[2, 0:2]); style(a)
    for i, (sa, c) in enumerate(zip(J['stage_a'], (C1, C2))):
        a.semilogy([o['d'] for o in sa['scan']], [o['mse'] for o in sa['scan']], color=c, lw=1.6, label='%s (minimum at %.0f nm)' % (sa['sheet'], sa['d']))
    a.axvspan(100, 250, color=MUTED, alpha=0.12, lw=0); a.text(175, 2.2, 'expected range', ha='center', fontsize=7.5, color=INK2)
    a.set_xlabel('imposed ZnS thickness (nm)', **lab); a.set_ylabel('MSE, 430-1000 nm', **lab); a.set_ylim(1.5, 1500); a.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='lower right')
    a.set_title('Thickness scan in the transparent range (Sellmeier index refitted at each thickness)', **ttl)
    a = fig.add_subplot(gs[2, 2]); style(a); a.fill_between(wl, Nt.real, Nb.real, color=C1, alpha=0.2, lw=0, label='top to bottom of the film'); a.plot(wl, N.real, color=C1, lw=1.8, label='n, depth average')
    a.set_xlabel('wavelength (nm)', **lab); a.set_ylabel('n', **lab); a.legend(frameon=False, fontsize=7.5, labelcolor=INK2); a.set_title('Refractive index', **ttl)
    a = fig.add_subplot(gs[2, 3]); style(a); a.semilogy(wl, np.maximum(N.imag, 1e-5), color=C1, lw=1.8); a.set_ylim(1e-4, 3); a.set_xlim(205, 600)
    a.set_xlabel('wavelength (nm)', **lab); a.set_ylabel('k', **lab); a.set_title('Extinction coefficient', **ttl)
    fig.savefig(OUT + 'zns_summary.png', dpi=180, facecolor=SURF, bbox_inches='tight')
    print('->', OUT)
    for r in (base, gr):
        print('+'.join(r['opts']) or 'homogeneous', 'MSE', [round(m, 2) for m in r['mse']], {k: '%.2f +- %.2f' % (r['struct'][k], r['ci'][k]) for k in r['names']})
    print('n,k (graded, depth average):', {w: (round(a_, 3), round(b_, 4)) for w, (a_, b_) in gr['nk'].items()})


if __name__ == '__main__':
    main()
