#!/usr/bin/env python3
"""Compare the author's LightTools BSDF of the hexagonal hemisphere MLA (lt_hemisphere_bsdf.mat,
BSDF_MLA[180, 90, 15], n_MLA = 1.30:0.05:2.00) with this repository's ray tracer, at n = 1.5
and 1.8: B_T(theta), p, first-pass escape of P_sub, and eta_ext through the same series."""
import numpy as np, scipy.io as sio, sys, os
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concurrent.futures import ProcessPoolExecutor
import cps2, materials as M, raytrace as rt, series

AL = 0.958336812040867 + 6.68678039868499j
TH = np.arange(90) + 0.5; W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
B = sio.loadmat('lt_hemisphere_bsdf.mat')['BSDF_MLA']
NS = np.round(np.arange(1.30, 2.001, 0.05), 2)


def lt_bsdf(n):
    k = int(np.argmin(abs(NS - n)))
    b = B[:, :, k]
    BT = b[:90, :].sum(0)                # transmitted rows 1:90 (MATLAB) -> escaped fraction per incidence bin
    BR = b[90:180, :][::-1, :]           # BSDF(180:-1:91, :) as the preprint script does: row 0 = 0-1 deg out
    return BT, BR, b


def dev(refl, n_sub):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=n_sub)


def _one_bin(a):
    n, i, th = a
    rng = np.random.default_rng([17, i, 21])
    esc, rc = rt._trace_chunk(1.0, n, 30000, 80, rng, float(th), 'cap')
    good = np.isfinite(rc); ang = np.degrees(np.arccos(np.clip(rc[good], 0.0, 1.0)))
    h, _ = np.histogram(ang, bins=np.linspace(0.0, 90.0, 91)); return i, esc.mean(), h / 30000


def my_bsdf(n, ex):
    BT = np.zeros(90); BR = np.zeros((90, 90))
    for i, t, col in ex.map(_one_bin, [(n, i, TH[i], ) for i in range(90)]):
        BT[i] = t; BR[:, i] = col
    return BT, BR


def report(tag, BT, BR, n_sub):
    out = ['%-14s p(cos.sin) %.4f  closure(min/max col sum) %.4f/%.4f' % (tag, np.sum(BT * W) / np.sum(W), (BT + BR.sum(0)).min(), (BT + BR.sum(0)).max())]
    out.append('   B_T at 0.5/10/20/30/40/50/60/70/80 deg: ' + ' '.join('%.3f' % BT[i] for i in (0, 10, 20, 30, 40, 50, 60, 70, 80)))
    for name, refl in (('Al', AL), ('Ag', M.AG)):
        S = dev(refl, n_sub); r = cps2.solve_pol(S, npts=12000)
        es = r['air'] + r['sub']; R = cps2.stack_reflectance(S, TH); P = cps2.sub_angular(S, TH)
        v0 = P / P.sum(); e1 = BT @ v0; v1 = BR @ v0; mth = np.sum(TH * v1) / v1.sum()
        e, _ = series.eta_ext(BT, BR, R, P, n_term=80)
        A = 1 - np.sum(R * W) / np.sum(W); pc = np.sum(BT * W) / np.sum(W)
        out.append('   %s: eta_sub %.4f  A\' %.4f | 1st-pass escape of P_sub %.4f | <theta> of returned light %.1f deg | eta_ext series %.4f (closed form %.4f) | EQE %.4f'
                   % (name, es, A, e1, mth, e, pc / (pc + (1 - pc) * A), es * e))
    print('\n'.join(out), flush=True)


if __name__ == '__main__':
    b = B[:, :, 10]
    print('LightTools slice n=1.80: column 0 sums rows 0:90 = %.4f, rows 90:180 = %.4f; column 80: %.4f / %.4f' % (b[:90, 0].sum(), b[90:, 0].sum(), b[:90, 80].sum(), b[90:, 80].sum()))
    print('first rows of column 0 (transmitted part):', np.round(b[:6, 0], 4), '... reflected part rows 90..95:', np.round(b[90:96, 0], 4), 'rows 174..179:', np.round(b[174:180, 0], 4))
    Z = np.load('mla_bsdf.npz'); ar = np.asarray(Z['aspect'], float); k = int(np.argmin(abs(ar - 1.0)))
    with ProcessPoolExecutor(4) as ex:
        mine15 = my_bsdf(1.5, ex)
    for n_sub in (1.5, 1.8):
        print('\n===== n_sub = n_MLA = %.1f =====' % n_sub)
        BT, BR, _ = lt_bsdf(n_sub)
        report('LightTools', BT, BR, n_sub)
        if n_sub == 1.8:
            report('this tracer', Z['BT'][k], Z['BR'][k], n_sub)
        else:
            report('this tracer', mine15[0], mine15[1], n_sub)
    np.savez_compressed('lt_vs_tracer_hemisphere.npz', lt15=lt_bsdf(1.5)[0], lt18=lt_bsdf(1.8)[0], mine15=mine15[0], mine18=Z['BT'][k],
                        lt15_BR=lt_bsdf(1.5)[1], lt18_BR=lt_bsdf(1.8)[1], mine15_BR=mine15[1], mine18_BR=Z['BR'][k])
