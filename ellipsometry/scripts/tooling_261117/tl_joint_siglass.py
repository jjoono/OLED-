"""Si + glass joint fit, the way the lab's multi-sample snapshots do it
(NPB_tooling_S8_S13: Biaxial film / NTVE_JAW / EMA substrate, per-sample
thickness).  Reference numbers to check a CompleteEASE joint fit against.

Per material, every continuous-film piece on Si and on glass together:
  NPB    1_21 1_22 1_23 (Si)  1_31 1_32 1_33 (glass)
  TAPC   1_24 (Si)            1_34 (glass)
  Bphen  2_21 2_22 (Si)       2_31 2_32 (glass)
Shared: the film's Gen-Osc (the constrained organic rules of section 6: Einf 1,
UV pole 6-12 eV, Gaussians amp >= 0, Br 0.10-1.20 eV, centre above the onset,
sub-gap eps2 penalised), and for 'biaxial' the axis scales.
Per piece: thickness; on glass also the collected backside fraction f.
Si: NTVE_JAW fixed at 1.5 nm (the lab snapshot carries ~1.45 nm).
Glass: this set's glass (tl_ce_glass2 Cauchy), fixed; backside on.

  iso       one isotropic Gen-Osc
  biaxial   Gaussian amplitudes x a_x, a_y, a_z (mean fixed at 1, so the
            Gen-Osc keeps the overall band strength); UV pole isotropic.
            x = plane of incidence, y = along s, z = normal: all one azimuth
            can tell apart.  On glass the weak backside beams go through the
            axis-averaged film (the approximation of tl_bphen10)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_genosc as GO, tl_aniso as A, tl_bphen10 as B
import jnk_model as M, jnk_ref as R
from tmm_fix import tmm
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
HC = 1239.841984
DOX = 1.5
SETS = {'NPB': ['1_21', '1_22', '1_23', '1_31', '1_32', '1_33'],
        'TAPC': ['1_24', '1_34'],
        'Bphen': ['2_21', '2_22', '2_31', '2_32']}
is_glass = lambda sh: sh.split('_')[1].startswith('3')


def load(sh):
    g = is_glass(sh)
    wl, P, D, ang = L.se(sh, 300.0 if g else 245.0, 1000.0, 3 if g else 2)
    return wl, L.ncs(P, D), ang


def ncs_glass_biax(wl, ang, ex, ey, ez, d, ng, f, nb=20):
    air = R.n_air(wl) + 0j; a2 = air**2; g2 = ng**2
    Niso = B.sqrtc((ex + ey + ez) / 3)
    out = []
    for th in ang:
        rp = B.tmm_uni(wl, [a2, ex, g2], [a2, ez, g2], [d], th)[0]
        rs = B.tmm_uni(wl, [a2, ey, g2], [a2, ey, g2], [d], th)[1]
        Ip, Is, X = np.abs(rp)**2, np.abs(rs)**2, np.conj(rp) * rs
        if f > 0:
            thg = np.rad2deg(np.arcsin(np.sin(np.deg2rad(th)) / ng.real))
            tdp, tds = M.t_amp(wl, [air, Niso, ng], [d], th)
            tup, tus = M.t_amp(wl, [ng, Niso, air], [d], thg)
            rip, ris = tmm(wl, [ng, Niso, air], [d], thg)
            rbp, rbs = tmm(wl, [ng, air], [], thg)
            ap = tdp * rbp * tup; as_ = tds * rbs * tus
            for _ in range(nb):
                Ip = Ip + f * np.abs(ap)**2; Is = Is + f * np.abs(as_)**2; X = X + f * np.conj(ap) * as_
                ap = ap * rip * rbp; as_ = as_ * ris * rbs
        I = Ip + Is
        out.append(((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    N, C_, S_ = (np.stack([o[k] for o in out], 1) for k in range(3))
    P = np.sqrt(N**2 + C_**2 + S_**2)               # 'normalised' depolarization convention (section 4)
    return N / P, C_ / P, S_ / P


def eps_axes(p, wl, a):
    out = []
    for ai in a:
        q = np.array(p, float); q[2::3] *= ai; out.append(GO.eps_genosc(q, wl))
    return out


def fit(mat, mode):
    shs = SETS[mat]; data = [load(sh) for sh in shs]
    gl = json.load(open(W + 'ce_glass2.json')); Ac, Bc = gl['cauchy']
    p0 = np.array(json.load(open(W + 'ce_genosc.json'))[mat]['pick']['p'], float)
    ng_osc = (len(p0) - 2) // 3; lo_p, hi_p = GO.bounds(mat, ng_osc); npar = len(p0)
    na = 2 if mode == 'biaxial' else 0
    glass_idx = [j for j, sh in enumerate(shs) if is_glass(sh)]
    d0 = [json.load(open(W + 'ce_glass2.json'))['d'].get(sh) or json.load(open(W + 'ce_thickness_ntve0.json'))[sh]['d_ntve0'] for sh in shs]
    f0 = [gl['f'][shs[j]] for j in glass_idx]
    x0 = np.r_[p0, [1.0] * na, d0, f0]
    LO = np.r_[lo_p, [0.3] * na, [1.0] * len(shs), [0.0] * len(glass_idx)]
    HI = np.r_[hi_p, [1.45] * na, [200.0] * len(shs), [1.0] * len(glass_idx)]
    eon = HC / GO.ONSET_NM[mat]; Eg = np.linspace(1.24, eon - 0.10, 40)

    def piece(x, j):
        p = x[:npar]; a = (x[npar], x[npar + 1], 3 - x[npar] - x[npar + 1]) if na else (1.0, 1.0, 1.0)
        d = x[npar + na + j]; wl, meas, ang = data[j]
        ex, ey, ez = eps_axes(p, wl, a)
        if is_glass(shs[j]):
            f = x[npar + na + len(shs) + glass_idx.index(j)]
            ng = Ac + Bc / (wl / 1000)**2 + 0j
            N, C_, S_ = ncs_glass_biax(wl, ang, ex, ey, ez, d, ng, f)
        else:
            air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
            N, C_, S_ = A.ncs_biaxial(wl, ang, [air, ex, ox, si], [air, ey, ox, si], [air, ez, ox, si], [d, DOX])
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])

    def res(x):
        return np.concatenate([piece(x, j) for j in range(len(shs))] + [C.PEN_W * GO.eps_genosc(x[:npar], HC / Eg).imag])

    r = least_squares(res, np.clip(x0, LO + 1e-6, HI - 1e-6), bounds=(LO, HI), x_scale='jac', max_nfev=600)
    x = r.x; J = r.jac; s2 = np.sum(r.fun**2) / max(r.fun.size - x.size, 1)
    ci = 1.645 * np.sqrt(np.clip(np.diag(np.linalg.pinv(J.T @ J) * s2), 0, None))
    per = {sh: float(M.mse_ncs(piece(x, j), (npar + na) / len(shs) + 1 + (1 if is_glass(sh) else 0))) for j, sh in enumerate(shs)}
    wl = np.array([350.0, 450.0, 550.0, 633.0])
    a = (x[npar], x[npar + 1], 3 - x[npar] - x[npar + 1]) if na else (1.0, 1.0, 1.0)
    nk = {}
    for lab, e in zip(('x', 'y', 'z'), eps_axes(x[:npar], wl, a)):
        n = B.sqrtc(e); nk[lab] = {str(int(w)): [float(v.real), float(v.imag)] for w, v in zip(wl, n)}
    return dict(mat=mat, mode=mode, pieces=shs, p=[float(v) for v in x[:npar]], a=[float(v) for v in a],
                d={sh: float(x[npar + na + j]) for j, sh in enumerate(shs)},
                d_ci90={sh: float(ci[npar + na + j]) for j, sh in enumerate(shs)},
                f={shs[j]: float(x[npar + na + len(shs) + k]) for k, j in enumerate(glass_idx)},
                mse=per, mse_mean=float(np.mean(list(per.values()))), nk=nk)


def main():
    from multiprocessing import Pool
    jobs = [(m, mode) for m in SETS for mode in ('iso', 'biaxial')]
    with Pool(4) as pool:
        res = pool.starmap(fit, jobs)
    json.dump(res, open(W + 'joint_siglass.json', 'w'), indent=1)
    for r in res:
        print('\n== %s, %s  (NTVE %.1f nm on Si; mean MSE %.2f)%s' % (r['mat'], r['mode'], DOX, r['mse_mean'],
              '   a_x / a_y / a_z = %.2f / %.2f / %.2f' % tuple(r['a']) if r['mode'] == 'biaxial' else ''))
        for sh in r['pieces']:
            print('   %s %-5s d = %6.2f +- %.2f nm   %s MSE %6.2f' % (sh, 'glass' if is_glass(sh) else 'Si', r['d'][sh], r['d_ci90'][sh],
                  ('f = %.2f ' % r['f'][sh]) if sh in r['f'] else '        ', r['mse'][sh]))
        print('   n,k (x axis): ' + '  '.join('%s: %.3f+%.3fj' % (w, *v) for w, v in r['nk']['x'].items()))
        if r['mode'] == 'biaxial':
            print('   n,k (z axis): ' + '  '.join('%s: %.3f+%.3fj' % (w, *v) for w, v in r['nk']['z'].items()))


if __name__ == '__main__':
    main()
