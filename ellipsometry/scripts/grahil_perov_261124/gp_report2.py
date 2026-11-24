"""Second figure of the GraHIL / perovskite session -> $ELLIPS_OUT/grahil_perov/grahil_perov_device.png

  a  angular emission: measured, the paper's optics, the two-layer structure with the scattered light
     followed (nominal geometry and the best geometry the micrograph allows), and its coherent part alone
  b  where the photons go, per description (direct air / trapped in glass / re-absorbed by perovskite / lost)
  c  layer2 fit quality with the surface roughness imposed (effective-medium layer; long-scale height spread)
  d  band edge of the dense perovskite: ellipsometry, luminescence-stitched, published FA/Cs lead bromide
     shifted to this film's edge, and the film-average k of the paper's absorbance"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, csv, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_device_check as DC
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

W = _os.path.join(ELLIPS_OUT, ''); OUT = _os.path.join(ELLIPS_OUT, 'grahil_perov', '')
INK, INK2, MUTED, GRID, AXIS, SURF = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
C1, C2, C3, C4 = '#2a78d6', '#eb6834', '#1baf7a', '#b05cc5'          # categorical slots 1-4
DASGUPTA = _os.environ.get('DASGUPTA_CSV', '')


def style(a):
    a.set_facecolor(SURF); a.grid(color=GRID, lw=0.6); a.tick_params(colors=MUTED, labelsize=8)
    for s in a.spines.values(): s.set_color(AXIS)


def pick(rows, sub='lost', scatter='follow', sig=15.0, mode='volume'):
    return [r for r in rows if r['sub'] == sub and r['scatter'] == scatter and r['sig'] == sig and r['x_d_mode'] == mode][0]


def main():
    M = json.load(open(W + 'device_model.json')); R = M['results']; th = np.array(M['theta']); Im = np.array(M['measured'])
    F = json.load(open(W + 'device_fit.json'))['families'] if _os.path.exists(W + 'device_fit.json') else None
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 9.2), facecolor=SURF); fig.subplots_adjust(hspace=0.36, wspace=0.26, left=0.07, right=0.98, top=0.93, bottom=0.07)
    lab = dict(color=INK2, fontsize=8.5); ttl = dict(color=INK, fontsize=9.5, loc='left')
    # a
    a = ax[0, 0]; style(a); fine = np.linspace(0, 85, 86); a.plot(fine, np.cos(np.deg2rad(fine)), color=MUTED, lw=1.1, label='Lambertian')
    p0 = pick(R['one layer, paper optics']); a.plot(th, p0['I'], color=C2, lw=1.8, marker='s', ms=4, mec=SURF, mew=0.8, label='one layer, n = 2.3 (the paper\'s optics), rms %.3f' % p0['rms'])
    p1 = pick(R['two layers']); a.plot(th, p1['I'], color=C1, lw=1.8, marker='o', ms=4, mec=SURF, mew=0.8, label='two layers, scattered light followed, rms %.3f' % p1['rms'])
    if F:
        b = F['two, followed'][0]['best']; a.plot(th, b['I'], color=C3, lw=1.8, marker='^', ms=4.5, mec=SURF, mew=0.8, label='the same, best geometry within the micrograph, rms %.3f' % b['rms'])
    p2 = pick(R['two layers'], scatter='lost'); a.plot(th, p2['I'], color=C1, lw=1.2, ls=(0, (4, 2)), label='two layers, coherent part alone, rms %.3f' % p2['rms'])
    a.plot(th, Im, color=INK, lw=0, marker='o', ms=7, mfc='none', mew=1.6, label='measured (Extended Data Fig. 7g)')
    a.set_xlabel('emission angle (deg)', **lab); a.set_ylabel('EL intensity / normal', **lab); a.set_ylim(0, 1.2); a.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='lower left')
    a.set_title('Angular emission of the green PeLED', **ttl)
    # b
    a = ax[0, 1]; style(a); a.grid(axis='x', visible=False)
    rows = [('one layer\n(paper optics)', pick(R['one layer, paper optics'])), ('two layers\nno scattering', pick(R['two layers, no scattering'])),
            ('two layers\nscatter = loss', pick(R['two layers'], scatter='lost')), ('two layers\nscatter followed', pick(R['two layers'])),
            ('same, micrograph\nthicknesses', pick(R['SEM geometry, FTO undercoat']))]
    keys = [('air', 'into air', C1), ('sub', 'trapped in the glass', C3), ('Ap', 're-absorbed by the perovskite', C2), ('Apar', 'lost in FTO, HIL, Al', MUTED), ('S', 'scattered, not followed', C4)]
    x = np.arange(len(rows)); bot = np.zeros(len(rows))
    for k, name, col in keys:
        v = np.array([r['fates'][k] for _, r in rows])
        if v.max() < 1e-6: continue
        a.bar(x, v, 0.62, bottom=bot, color=col, edgecolor=SURF, lw=1.5, label=name); bot += v
    for i, (_, r) in enumerate(rows):
        a.text(i, 1.02, 'EQE at\n$\\eta_{rad}$=1: %.0f %%' % (100 * r['eqe_pl1']), ha='center', va='bottom', fontsize=7.5, color=INK2)
    a.set_xticks(x); a.set_xticklabels([n for n, _ in rows], fontsize=7.5, color=INK2); a.set_ylim(0, 1.22); a.set_ylabel('share of emitted photons', **lab)
    a.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3); a.set_title('Where the photons go (substrate light counted as lost)', **ttl)
    # c
    a = ax[1, 0]; style(a); J = json.load(open(W + 'rough_scan.json')); sc = {r['variant']: r for r in J['scan']}; b0 = J['base']
    dr = [(b0['struct']['d_r'], b0['mse_ncs'])] + [(sc[k]['struct']['d_r'], sc[k]['mse_ncs']) for k in ('d_r 3', 'd_r 6', 'd_r 10', 'd_r 15', 'd_r 20', 'd_r 30') if k in sc]
    sg = [(0.0, b0['mse_ncs'])] + [(sc[k]['struct']['sig_s'], sc[k]['mse_ncs']) for k in ('sig_s 5', 'sig_s 10', 'sig_s 13.5', 'sig_s 20') if k in sc]
    a.plot(*zip(*dr), color=C2, lw=1.8, marker='s', ms=4.5, mec=SURF, mew=0.8, label='effective-medium roughness layer (lateral scale << wavelength)')
    a.plot(*zip(*sg), color=C1, lw=1.8, marker='o', ms=4.5, mec=SURF, mew=0.8, label='height spread, long lateral scale (coherent average)')
    a.axvspan(3, 9, color=C2, alpha=0.10, lw=0); a.axvspan(10, 13.5, color=C1, alpha=0.10, lw=0)
    a.text(6, 46, 'micro-roughness\nfrom top-view SEM', ha='center', fontsize=7, color=C2); a.text(11.75, 34, 'undulation\nfrom the\ncross-section', ha='center', fontsize=7, color=C1)
    a.set_xlabel('imposed roughness (nm): layer thickness, or rms height', **lab); a.set_ylabel('MSE of the layer2 fit', **lab); a.set_ylim(0, 60); a.set_xlim(-0.5, 21)
    a.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='upper left'); a.set_title('Ellipsometry fit with the roughness imposed', **ttl)
    # d
    a = ax[1, 1]; style(a); st = np.loadtxt(OUT + 'Perov_dense_nk_stitched.csv', delimiter=',', skiprows=1); m = (st[:, 0] >= 500) & (st[:, 0] <= 600)
    a.semilogy(st[m, 0], np.maximum(st[m, 4], 1e-7), color=C1, lw=1.3, ls=(0, (4, 2)), label='ellipsometry fit (edge shape set by the basis)')
    a.semilogy(st[m, 0], st[m, 2], color=C1, lw=2.0, label='stitched: luminescence edge + Urbach tail')
    S = json.load(open(OUT + 'stitch.json'))
    if DASGUPTA and _os.path.exists(DASGUPTA) and 'published' in S:
        rows_ = list(csv.reader(open(DASGUPTA, newline=''))); h = rows_[0]; d = np.array([[float(v) if v not in ('', 'nan') else np.nan for v in rr] for rr in rows_[1:]])
        wd, kd = d[:, 0], d[:, h.index('k_Br_100%_FaCsPbBrI')]; sh = S['published']['edge_shift_nm']; mm = (wd - sh >= 500) & (wd - sh <= 600)
        a.semilogy(wd[mm] - sh, kd[mm], color=C3, lw=1.8, label='published Cs$_{0.17}$FA$_{0.83}$PbBr$_3$ (PDS-stitched), shifted %.0f nm' % (-sh))
    a.semilogy(DC.K_PAPER[0], np.maximum(DC.K_PAPER[1], 1e-7), color=C2, lw=0, marker='o', ms=4.5, mec=SURF, mew=0.8, label='film average from absorbance (paper SI Fig. 6a)')
    a.set_ylim(1e-6, 1.0); a.set_xlim(500, 600); a.set_xlabel('wavelength (nm)', **lab); a.set_ylabel('k', **lab); a.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='lower left')
    a.set_title('Band edge of the dense perovskite', **ttl)
    fig.savefig(OUT + 'grahil_perov_device.png', dpi=200, facecolor=SURF, bbox_inches='tight'); print('->', OUT + 'grahil_perov_device.png')


if __name__ == '__main__':
    main()
