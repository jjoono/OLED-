"""References without an angular filter.
(1) Hemispherical MLA (LightTools BSDF, n_MLA = n_sub), angle-preserving device: trans-scale series.
(2) Strong-microcavity bottom emitter: thin Ag anode (20 nm) instead of ITO, planar glass exit,
    ETL/HTL scanned; Al or Ag cathode.  Stored per wavelength so that the best cavity can be
    chosen for every spectrum and cone afterwards."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, scipy.io as sio, time, sys
from multiprocessing import Pool
import sel
from sel import cps2

BSDF = sys.argv[1] if len(sys.argv) > 1 else os.path.join(sel.HERE, 'data', 'lt_hemisphere_bsdf.mat')
DV = np.load('data/devices.npz')
EA = np.arange(91.0)                       # air bins 1 deg


def mla_ref():
    B = sio.loadmat(BSDF)['BSDF_MLA']      # [180, 90, n], n_MLA = 1.30:0.05:2.00
    out = {}
    for k in sel.STACKS:
        ns = float(DV[f'{k}_n_sub']); kn = int(round((ns - 1.30) / 0.05))
        b = B[:, :, kn]; BTm = b[:90]; BR = b[90:180][::-1]          # BTm[air, in], BR[sub_out, in]
        P = DV[f'{k}_P'].reshape(sel.LAM.size, 90, 2).sum(2)
        R = DV[f'{k}_R'].reshape(sel.LAM.size, 90, 2).mean(2)
        cl = np.zeros((sel.LAM.size, len(sel.BETAS)))
        for i in range(sel.LAM.size):
            x = np.linalg.solve(np.eye(90) - R[i][:, None] * BR, P[i])
            air = BTm @ x
            for j, bt in enumerate(sel.BETAS):
                f = np.clip(bt - EA[:-1], 0, 1)
                cl[i, j] = air @ f
        out[k] = cl
    return out


def cav_one(args):
    cath, d_etl, d_htl = args
    es, cl = [], []
    for j, lam in zip(sel.IL, sel.LAM):
        nk = lambda m: sel.nk(m, j)
        org = lambda m, d: (nk(m + '_o').real, nk(m + '_e').real, d)
        S = cps2.Stack(float(lam), (nk('EML_o').real, nk('EML_e').real), 25.0, 12.5,
                       above=[org('B3', d_etl), (nk(cath), nk(cath), 100.0)],
                       below=[org('TCTA', 10.0), org('TAPC', d_htl), (nk('Ag'), nk('Ag'), 20.0)], n_sub=1.51, h=0.865)
        r = cps2.solve_pol(S, npts=3000)
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(sel.NB, 10).sum(1)
        R = np.clip(np.real(cps2.stack_reflectance(S, sel.TH)), 0, 1)
        dev = dict(P=(pf / pf.sum())[None], R=R[None])
        Rf, Tf = sel.tmm_RT([], float(lam), sel.TH, 1.51)
        air = Tf * dev['P'][0] / (1 - Rf * R)                       # planar, angle preserving
        es.append(r['air'] + r['sub']); cl.append(air @ sel.cone_frac(1.51).T)
    return args, np.array(es), np.array(cl)


if __name__ == '__main__':
    t0 = time.time()
    m = mla_ref(); np.savez('data/ref_mla.npz', **m); print('MLA done', time.time() - t0, flush=True)
    E = np.arange(20.0, 301.0, 20.0); H = np.arange(20.0, 301.0, 20.0)
    jobs = [(c, e, h) for c in ('Ag', 'Al') for e in E for h in H]
    res = {}
    with Pool(4) as p:
        for c, (a, es, cl) in enumerate(p.imap_unordered(cav_one, jobs, chunksize=2)):
            res[a] = (es, cl)
            if c % 40 == 0:
                print(c, len(jobs), f'[{time.time() - t0:.0f}s]', flush=True)
    keys = sorted(res)
    np.savez('data/ref_cavity.npz', keys=np.array([f'{c}|{e}|{h}' for c, e, h in keys]),
             eta_sub=np.array([res[k][0] for k in keys]), cones=np.array([res[k][1] for k in keys]))
    print('done', time.time() - t0)
