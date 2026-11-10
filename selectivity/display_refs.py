"""Display metrics (on-axis cd/A, luminance and u'v' vs angle, EQE) for the non-cavity configurations."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, numpy as np, scipy.io as sio
import sel, common, optfilter as O, run_display as RD
from sel import cps2

DV = RD.DV; SL = RD.SL
B = sio.loadmat('data/lt_hemisphere_bsdf.mat')['BSDF_MLA']


def dev(st):
    return dict(eta_sub=DV[f'{st}_eta_sub'], P=DV[f'{st}_P'], R=DV[f'{st}_R'], n_sub=float(DV[f'{st}_n_sub']))


def mla_I(st):
    d = dev(st); kn = int(round((d['n_sub'] - 1.30) / 0.05)); b = B[:, :, kn]; BTm = b[:90]; BR = b[90:180][::-1]
    P = d['P'].reshape(-1, 90, 2).sum(2); R = d['R'].reshape(-1, 90, 2).mean(2)
    air = np.array([BTm @ np.linalg.solve(np.eye(90) - R[i][:, None] * BR, P[i]) for i in range(sel.LAM.size)])
    return RD.rebin(air, np.arange(90.0), np.arange(1.0, 91.0)) / RD.SR[None, :]


def filt_I(st, f, rnd):
    d = dev(st); Rf, Tf = sel.filter_RT(f, d['n_sub'])
    return RD.to_air_angle(sel.solve(d, Rf, Tf, *rnd), d['n_sub'])


def thin_ag_cavity(cath, e, h):
    es, P, R = [], [], []
    for lam in sel.LAM:
        j = int(lam) - 400; nk = lambda m: sel.nk(m, j); org = lambda m, d: (nk(m + '_o').real, nk(m + '_e').real, d)
        S = cps2.Stack(float(lam), (nk('EML_o').real, nk('EML_e').real), 25.0, 12.5, above=[org('B3', e), (nk(cath), nk(cath), 100.0)],
                       below=[org('TCTA', 10.0), org('TAPC', h), (nk('Ag'), nk('Ag'), 20.0)], n_sub=1.51, h=0.865)
        r = cps2.solve_pol(S, npts=2000)
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(sel.NB, 10).sum(1)
        es.append(r['air'] + r['sub']); P.append(pf / pf.sum()); R.append(np.clip(np.real(cps2.stack_reflectance(S, sel.TH)), 0, 1))
    return dict(eta_sub=np.array(es), P=np.array(P), R=np.array(R), n_sub=1.51)


if __name__ == '__main__':
    CAV = np.load('data/ref_cavity.npz'); out = {'angles': RD.AC.tolist(), 'configs': {}}
    slab = tuple(SL['g0.8_t1.0_w1'])
    for sp in RD.SPECS:
        spec = common.SPECTRA[sp]
        cfg = {}
        cfg['Al planar (ITO 150)'] = (RD.to_air_angle(sel.solve(dev('Al'), *sel.filter_RT('fresnel', 1.51), *sel.mixer(0.0)), 1.51), dev('Al')['eta_sub'])
        cfg['low-loss n1.8 + MLA'] = (mla_I('LL18'), dev('LL18')['eta_sub'])
        E = np.einsum('l,kl,klb->kb', spec, CAV['eta_sub'], CAV['cones'])
        ks = [i for i, k in enumerate(CAV['keys']) if k.startswith('Ag|')]; kb = CAV['keys'][max(ks, key=lambda i: E[i, 0])]
        _, e, h = kb.split('|'); d = thin_ag_cavity('Ag', float(e), float(h))
        cfg[f'thin-Ag strong cavity (ETL {float(e):.0f}, HTL {float(h):.0f})'] = (RD.to_air_angle(RD.planar_air(d), 1.51), d['eta_sub'])
        best = None
        for mat in ('TiO2/SiO2', 'ZnS/LiF'):
            P = O.Problem('LL18', sp, 5.0, mat, key='fwd')
            for N in (8, 12, 16):
                for lb in np.arange(540., 760., 5.):
                    v = P.eqe(O.edge_seq(mat, N, lb))
                    if best is None or v > best[0]: best = (v, mat, N, float(lb))
        cfg[f'low-loss n1.8 + scatterer + edge DBR ({best[1]} {best[2]}p, {best[3]:.0f} nm)'] = (filt_I('LL18', ('edge',) + best[1:], slab), dev('LL18')['eta_sub'])
        for b in (20.0, 10.0):
            cfg[f'low-loss n1.8 + scatterer + ideal {b:g} deg filter'] = (filt_I('LL18', ('ideal', b), slab), dev('LL18')['eta_sub'])
        for k, (I, es) in cfg.items():
            m = RD.display_metrics(I, es, spec); out['configs'][f'{sp}|{k}'] = m
            print(f"{sp:13s} {k:62s} cd/A {m['cdA0']:6.0f}  EQE {m['EQE']:.3f}  half-lum {m['half_lum_angle']:4.1f} deg  du'v'(<=30) {m['duv30']:.4f}  (<=60) {m['duv60']:.4f}", flush=True)
    json.dump(out, open('data/display_refs.json', 'w'))
