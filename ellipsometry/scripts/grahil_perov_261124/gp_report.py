"""Deliverables of the GraHIL / perovskite-NC analysis -> $ELLIPS_OUT/grahil_perov/

  Perov_dense_nk.csv     dense perovskite n,k (base fit) on a 1 nm grid, with the envelope over
                         all fit variants, and the effective n,k of the lower layer
  GraHIL_nk.csv          GraHIL n,k (written by gp_grahil.py, copied here)
  lib_*.csv              the same media as wavelength_nm,n,k (the layout of the n,k library)
  structure.json         layer thicknesses and the numbers quoted in the README
  grahil_perov_summary.png   fit, n,k, depth profile, device angular emission"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, shutil, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_load as G, gp_model as P, gp_perov as PV, gp_device_check as DC
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

W = _os.path.join(ELLIPS_OUT, ''); OUT = _os.path.join(ELLIPS_OUT, 'grahil_perov', '')
INK, INK2, MUTED, GRID, AXIS, SURF = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
C1, C2, C3 = '#2a78d6', '#eb6834', '#1baf7a'          # validated categorical slots 1-3 (blue, orange, aqua)
SC_HOLD = 450.0                                       # nm; below this the scattering extinction of the lower layer is held constant


def style(a):
    a.set_facecolor(SURF); a.grid(color=GRID, lw=0.6); a.tick_params(colors=MUTED, labelsize=8)
    for s in a.spines.values(): s.set_color(AXIS)


def main():
    _os.makedirs(OUT, exist_ok=True)
    J = json.load(open(W + 'perov_fits.json')); base = [r for r in J['variants'] if r['variant'] == 'base'][0]
    wl = np.arange(205.0, 1001.0, 1.0)
    N0 = PV.nk_of(base, wl); allN = np.array([PV.nk_of(r, wl) for r in J['variants']])
    v = base['struct']; one = np.ones(len(wl), complex)
    # the lambda^-p law is seen only where light reaches the lower layer (the dense layer is opaque
    # below ~450 nm) and cannot keep growing into the UV: held at its 450 nm value below that
    Nlow = P.bruggeman(N0, one, v['f_b']) + 1j * v['k_sc'] * (633.0 / np.maximum(wl, SC_HOLD))**v['p']
    cols = np.c_[wl, N0.real, N0.imag, allN.real.min(0), allN.real.max(0), allN.imag.min(0), allN.imag.max(0), Nlow.real, Nlow.imag]
    np.savetxt(OUT + 'Perov_dense_nk.csv', cols, delimiter=',', fmt='%.5f', comments='',
               header='wavelength_nm,n_dense,k_dense,n_min,n_max,k_min,k_max,n_lower_layer,k_lower_layer_incl_scattering')
    if _os.path.exists(W + 'GraHIL_nk.csv'): shutil.copy(W + 'GraHIL_nk.csv', OUT + 'GraHIL_nk.csv')
    # the same three media in the n,k-library layout (wavelength_nm,n,k) for TMM / dipole codes
    lib = lambda name, N: np.savetxt(OUT + name, np.c_[wl, N.real, N.imag], delimiter=',', fmt='%.5f', comments='', header='wavelength_nm,n,k')
    lib('lib_Perov_dense.csv', N0); lib('lib_Perov_lower_with_scattering.csv', Nlow); lib('lib_Perov_lower_no_scattering.csv', P.bruggeman(N0, one, v['f_b']))
    lib('lib_GraHIL.csv', PV.grahil_N(wl)[0])
    gf = json.load(open(W + 'grahil_fits.json')); gi = [r for r in gf['variants'] if r['variant'] == 'iso'][0]
    dev = json.load(open(W + 'device_check.json'))
    q = [450, 500, 520, 530, 540, 550, 600, 633, 800, 1000]; iq = [int(w - 205) for w in q]
    summ = dict(
        grahil=dict(thickness_nm=gi['d'], n_550=gi['n_o']['550'][0], k_633=gi['n_o']['633'][1], mse=gi['mse'],
                    thickness_range_oxide_0p5_to_3nm=[min(r['d'] for r in gf['oxide_scan']), max(r['d'] for r in gf['oxide_scan'])]),
        perovskite_on_si=dict(model='EMA roughness / dense / lower (perovskite + void, scattering) / GraHIL %.1f nm / SiO2 %.1f nm / Si' % (PV.D_G, PV.D_OX),
                              base=dict(struct=base['struct'], mse_ncs=base['mse_ncs'], bands=base['bands']),
                              variants={r['variant']: dict(mse_ncs=r['mse_ncs'], rms_R=r['rms_R'], struct=r['struct']) for r in J['variants']},
                              total_perovskite_nm=v['d_t'] + v['d_b']),
        dense_nk={str(w): dict(n=float(N0[i].real), k=float(N0[i].imag), n_range=[float(allN.real[:, i].min()), float(allN.real[:, i].max())],
                               k_range=[float(allN.imag[:, i].min()), float(allN.imag[:, i].max())]) for w, i in zip(q, iq)},
        lower_layer_nk={str(w): [float(Nlow[i].real), float(Nlow[i].imag)] for w, i in zip(q, iq)},
        lower_layer_note='effective medium of the base fit: Bruggeman(dense, void, f_b) + i k_sc (633/lambda)^p; the extinction is scattering '
                         '(light leaves the specular beam), determined over ~480-1000 nm and held constant below %g nm' % SC_HOLD,
        device_angular={r['label'].strip(): dict(rms=r['rms'], I=r['I']) for r in dev['rows']}, device_theta=dev['theta'], device_measured=dev['measured'],
        device_reference=dict(lambertian=dev['lambertian'], fully_scattered=dev['diffuse']))
    json.dump(summ, open(OUT + 'structure.json', 'w'), indent=1)

    # ---------------------------------------------------------------- figure
    S = PV.setup(205.0, 1000.0, 1); Np = PV.nk_of(base, S['wl']); rp, rs = PV.stack_r(S, v, Np)
    ps, dl = P.psidelta(*P.ncs_from_r(rp, rs)); d = G.load('layer2', 205.0, 1000.0, 1)
    fig = plt.figure(figsize=(13.5, 9.6), facecolor=SURF); gs = fig.add_gridspec(3, 6, height_ratios=[1, 1, 1.25], hspace=0.42, wspace=0.9)
    for a in range(3):
        ax = fig.add_subplot(gs[0, 2 * a:2 * a + 2]); style(ax)
        ax.plot(S['wl'], d['psi'][:, a], color=INK2, lw=2.6, alpha=0.45, label='measured'); ax.plot(S['wl'], ps[:, a], color=C1, lw=1.3, label='model')
        ax.set_title('%g°' % S['ang'][a], color=INK, fontsize=9, loc='left'); ax.set_ylabel('Ψ (deg)', color=INK2, fontsize=8.5)
        if a == 0: ax.legend(frameon=False, fontsize=8, labelcolor=INK2)
        ax2 = fig.add_subplot(gs[1, 2 * a:2 * a + 2]); style(ax2)
        dm = d['del'][:, a] % 360; mm = dl[:, a].copy()
        for arr in (dm, mm): arr[np.r_[False, np.abs(np.diff(arr)) > 180]] = np.nan          # do not draw the 0/360 wraps
        ax2.plot(S['wl'], dm, color=INK2, lw=2.6, alpha=0.45); ax2.plot(S['wl'], mm, color=C1, lw=1.3)
        ax2.set_ylabel('Δ (deg)', color=INK2, fontsize=8.5); ax2.set_xlabel('wavelength (nm)', color=INK2, fontsize=8.5)
    # n, k
    ax = fig.add_subplot(gs[2, 0:2]); style(ax); m = (wl >= 400) & (wl <= 800)
    ax.fill_between(wl[m], allN.real.min(0)[m], allN.real.max(0)[m], color=C1, alpha=0.25, lw=0)
    ax.plot(wl[m], N0.real[m], color=C1, lw=2, label='n, dense perovskite (this work)')
    ax.plot([505, 575], [2.3, 2.3], color=C2, lw=2, label='n = 2.3 assumed in the paper')
    ax.plot(wl[m], Nlow.real[m], color=C3, lw=2, label='n, lower layer (effective)')
    ax.set_ylim(1.3, 2.5); ax.set_xlabel('wavelength (nm)', color=INK2, fontsize=8.5); ax.set_ylabel('n', color=INK2, fontsize=8.5)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='center right', bbox_to_anchor=(1.0, 0.43)); ax.set_title('Refractive index', color=INK, fontsize=9, loc='left')
    ax = fig.add_subplot(gs[2, 2:4]); style(ax); m = (wl >= 480) & (wl <= 600)
    ax.fill_between(wl[m], allN.imag.min(0)[m], allN.imag.max(0)[m], color=C1, alpha=0.25, lw=0)
    ax.plot(wl[m], N0.imag[m], color=C1, lw=2, label='k, dense perovskite (this work)')
    fa = (v['d_t'] + v['f_b'] * v['d_b']) / (v['d_t'] + v['d_b'])
    ax.plot(wl[m], fa * N0.imag[m], color=C1, lw=1.3, ls=(0, (4, 2)), label='the same x %.2f (film average)' % fa)
    ax.plot(DC.K_PAPER[0], DC.K_PAPER[1], color=C2, lw=0, marker='o', ms=4.5, mec=SURF, mew=0.8, label='k from absorbance (paper SI Fig. 6a)')
    ax.set_xlabel('wavelength (nm)', color=INK2, fontsize=8.5); ax.set_ylabel('k', color=INK2, fontsize=8.5); ax.set_ylim(0, 0.31)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='upper right'); ax.set_title('Extinction coefficient at the band edge', color=INK, fontsize=9, loc='left')
    # device
    ax = fig.add_subplot(gs[2, 4:6]); style(ax); th = np.array(dev['theta']); fine = np.linspace(0, 85, 86)
    ax.plot(fine, np.cos(np.deg2rad(fine)), color=MUTED, lw=1.2, label='Lambertian')
    ax.plot(fine, DC.diffuse(fine), color=MUTED, lw=1.2, ls=(0, (4, 2)), label='fully scattered')
    pick = {r['label'].strip(): r for r in dev['rows']}
    se = pick['SE structure (base fit)      spread 15']; pa = pick['paper optics                 spread 15']
    ax.plot(th, pa['I'], color=C2, lw=2, marker='s', ms=4, mec=SURF, mew=0.8, label='paper optics, rms %.3f' % pa['rms'])
    ax.plot(th, se['I'], color=C1, lw=2, marker='o', ms=4, mec=SURF, mew=0.8, label='this work, rms %.3f' % se['rms'])
    ax.plot(th, dev['measured'], color=INK, lw=0, marker='o', ms=7, mfc='none', mew=1.6, label='measured (ED Fig. 7g)')
    ax.set_xlabel('emission angle (deg)', color=INK2, fontsize=8.5); ax.set_ylabel('EL intensity / normal', color=INK2, fontsize=8.5); ax.set_ylim(0, 1.2)
    ax.legend(frameon=False, fontsize=7.5, labelcolor=INK2, loc='lower left'); ax.set_title('Green PeLED, angular emission', color=INK, fontsize=9, loc='left')
    fig.suptitle('Perovskite NC film on GraHIL on Si: three-angle fit 205-1000 nm (MSE %.1f), dense %.0f nm on a low-index scattering %.0f nm, GraHIL %.0f nm'
                 % (base['mse_ncs'], v['d_t'], v['d_b'], PV.D_G), color=INK, fontsize=10, x=0.02, ha='left', y=0.975)
    fig.savefig(OUT + 'grahil_perov_summary.png', dpi=200, facecolor=SURF, bbox_inches='tight')
    print('->', OUT); print(json.dumps(summ['dense_nk'], indent=0)[:900])


if __name__ == '__main__':
    main()
