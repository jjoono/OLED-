"""Lateral spread (point-spread function) of the angle-selective recycling OLED (optimised DBR-4 / DBR-6 structures,
mu_a d 0.01, S 30 / S 1 back diffuser) at 530 nm, Monte Carlo over the recycling passes:
  emitter -> diffuser (angle ~ P_s) -> redistribution column of B (weight = reflectance) -> up through the stack at theta_o:
  exit with T_f(theta_o) (recorded: radius, air angle), else reflected (R_f) back down; lateral steps: diffuser
  d_s (tan th_i + tan th_o), upper stack 2 sum d_l tan th_l (layer angles by Snell, DBR fully traversed).
Reports r50 / r90 of the on-axis light (air angle <= 11 deg) and of all emitted light, for diffuser thickness d_s.
Compared with the substrate-MLA device (0.5 mm glass): r90 ~ 4-9 d from the earlier Monte Carlo = 2-4.5 mm."""
import json, numpy as np
import ce_opt as C
D = C.D; j = 130; lam = 530.0
def run(fn, N, mua, S, g, d_s, n_ph=400000, seed=1):
    s = json.load(open(fn))['structure']; x = [s[k] for k in ('d_ETL', 'd_HTL', 'lam_c', 'd_IZO_front', 'd_IZO_back', 'n_s')]
    P, Ld, Pair, (Ps, Tf, Rf, ns) = C.point(('B', j, tuple(x), N, 3000))
    B = C.BZ['S%g_g%g_a%g' % (S, g, mua)]
    jc = int(round(x[2] - 400))
    lay = [(C.IZO[j].real, x[4]), (D.B3o[j].real, x[0]), (D.EMo[j].real, 25.0), (D.TPo[j].real, x[1]), (D.HAT[j].real, 5.0), (C.IZO[j].real, x[3])]
    lay += [(C.ZNS[j].real, x[2] / 4 / C.ZNS[jc].real), (C.LIF[j].real, x[2] / 4 / C.LIF[jc].real)] * N
    TB = np.radians(C.TB)
    def up(th):           # lateral shift for one pass through the upper stack at medium angle th
        sn = ns * np.sin(th); tot = 0
        for n, d in lay:
            st = np.clip(sn / n, 0, 0.999); tot += d * st / np.sqrt(1 - st ** 2)
        return tot
    UP = np.array([up(t) for t in TB]) / 1000.0; TS = np.tan(np.minimum(TB, np.radians(85)))   # um
    rng = np.random.default_rng(seed); p0 = Ps / Ps.sum()
    ib = rng.choice(45, n_ph, p=p0); w = np.full(n_ph, Ps.sum()); xy = np.zeros((n_ph, 2)); alive = np.ones(n_ph, bool)
    rec_r, rec_w, rec_on = [], [], []
    Bc = B / np.maximum(B.sum(0), 1e-12); Bcum = np.cumsum(Bc, 0); refl = B.sum(0)
    for it in range(400):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        th_i = ib[a]; w[a] *= refl[th_i]
        ob = (Bcum[:, th_i].T < rng.random(a.size)[:, None]).sum(1).clip(0, 44)
        ph = rng.uniform(0, 2 * np.pi, (a.size, 2))
        step = d_s * (TS[th_i][:, None] * np.c_[np.cos(ph[:, 0]), np.sin(ph[:, 0])] + TS[ob][:, None] * np.c_[np.cos(ph[:, 1]), np.sin(ph[:, 1])])
        xy[a] += step
        ph2 = rng.uniform(0, 2 * np.pi, a.size); u = np.c_[np.cos(ph2), np.sin(ph2)]
        xy[a] += UP[ob][:, None] * u
        rec_r.append(np.hypot(*xy[a].T)); rec_w.append(w[a] * Tf[ob]); rec_on.append(ob < 3)
        w[a] *= Rf[ob]; xy[a] += UP[ob][:, None] * u; ib[a] = ob
        alive[a[w[a] < 1e-4]] = False
    r = np.concatenate(rec_r); wt = np.concatenate(rec_w); on = np.concatenate(rec_on)
    def q(m, f):
        o = np.argsort(r[m]); c = np.cumsum(wt[m][o]); return r[m][o][np.searchsorted(c, f * c[-1])]
    return q(on, 0.5), q(on, 0.9), q(np.ones_like(on), 0.5), q(np.ones_like(on), 0.9), float(np.sum(wt[on]) / np.sum(wt))
if __name__ == '__main__':
    for fn, N in (('results_ce/B_N4_mua0.01.json', 4), ('results_ce/B_N6_mua0.01.json', 6)):
        for S, g, d_s in ((30.0, 0.0, 1.0), (30.0, 0.0, 3.0), (1.0, 0.0, 0.5), (1.0, 0.0, 1.0)):
            r50, r90, a50, a90, fon = run(fn, N, 0.01, S, g, d_s)
            print('DBR%d  diffuser S %g g %g d_s %.1f um : on-axis r50 %.1f um r90 %.1f um | all light r50 %.1f r90 %.1f um'
                  % (N, S, g, d_s, r50, r90, a50, a90), flush=True)
