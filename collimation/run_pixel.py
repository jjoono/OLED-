"""Pixel-level study: scatterer thickness / mean free path / anisotropy, gap, pixel pitch and
reflective bank walls, for the scatterer + DBR collimator (FWHM 20 nm design, beta = 20 deg)."""
import json, sys, time, numpy as np
import model as Mo, pixel_mc as PM

LAM = np.arange(470.0, 601.0, 5.0)
R = json.load(open('results.json'))


def setup(fwhm=20.0, beta=20.0):
    row = next(r for r in R['rows'] if r['fwhm'] == fwhm and r['a_sc'] == 0.0 and r['beta'] == beta)
    nH, nL = {'TiO2/SiO2': (2.40, 1.46), 'Ta2O5/SiO2': (2.10, 1.46), 'SiNx/SiO2': (1.95, 1.46)}[row['mat']]
    ns, ds = Mo.dbr_stack(row['lam_b'], row['N'], nH, nL)
    W = R['weak']
    es, Pc, Rd, Td = [], [], [], []
    for lam in LAM:
        e, P, Rr = Mo.device(lam, W['d_exit'], W['d_refl'], cathode='ito')
        es.append(e); Pc.append(np.cumsum(P) / P.sum()); Rd.append(Rr); Td.append(Mo.tmm_T(ns, ds, lam, Mo.THR))
    sp = Mo.spectrum(LAM, 530.0, fwhm); es = np.array(es)
    wl = sp * es; eta_sub = wl.sum()
    return dict(row=row, lamcdf=np.cumsum(wl) / wl.sum(), Pcdf=np.array(Pc), Rdev=np.array(Rd), Tdbr=np.array(Td),
                eta_sub=eta_sub, ref=row['ref_eta_beta'], beta=beta)


def simulate(S, n, P, s, d, ls, g, alb, seed=1):
    t = s + d
    st, ang, hits, cx, cy, wt, lw1, lw2 = PM.run(n, seed, S['lamcdf'], S['Pcdf'], S['Rdev'], S['Tdbr'],
                                       P, s, t, ls if d > 0 else 0.0, g, alb, 20000)
    esc = st == 0
    inb = esc & (ang * 0.25 < S['beta'])
    own = (cx == 0) & (cy == 0)
    out = dict(lost=float(np.mean(st == 2)))
    # no walls: open plane
    out['open_total'] = S['eta_sub'] * float(np.sum(wt[esc]) / n)
    out['open_beta'] = S['eta_sub'] * float(np.sum(wt[inb]) / n)
    out['open_beta_own'] = S['eta_sub'] * float(np.sum(wt[inb & own]) / n)
    out['crosstalk'] = 1 - float(np.sum(wt[esc & own]) / max(np.sum(wt[esc]), 1e-12))
    out['mean_hits'] = float(np.mean(hits[esc])) if esc.any() else 0.0
    for Rw in (0.90, 0.95, 0.98, 0.995):
        f = Rw ** hits.astype(float)
        out[f'bank{Rw}_beta'] = S['eta_sub'] * float(np.sum((wt * f)[inb]) / n)
        out[f'bank{Rw}_total'] = S['eta_sub'] * float(np.sum((wt * f)[esc]) / n)
    for nm, lw in (('tir_air', lw1), ('tir_n125', lw2)):
        f = np.exp(lw)
        out[nm + '_beta'] = S['eta_sub'] * float(np.sum((wt * f)[inb]) / n)
        out[nm + '_total'] = S['eta_sub'] * float(np.sum((wt * f)[esc]) / n)
    return out


if __name__ == '__main__':
    t0 = time.time(); S = setup()
    print(f"design {S['row']['mat']} N={S['row']['N']} lam_b={S['row']['lam_b']}, eta_sub={S['eta_sub']:.3f}, "
          f"closed-form eta_beta={S['row']['eta_beta']:.3f}, strong-cavity ref={S['ref']:.3f}", flush=True)
    rows = []
    # (1) validation: huge pixel, thick strongly scattering slab
    v = simulate(S, 200000, 1e6, 0.0, 20.0, 0.5, 0.0, 1.0)
    print('validation (ideal randomizer limit): open eta_beta %.3f  total %.3f  lost %.4f' % (v['open_beta'], v['open_total'], v['lost']), flush=True)
    # (2) scatterer study, open plane (thickness d, mean free path ls, anisotropy g, albedo)
    for g in (0.0, 0.6, 0.8):
        for ls in (0.5, 1.0, 2.0):
            for d in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
                for alb in (1.0, 0.999, 0.995):
                    if d == 0.0 and (ls != 0.5 or alb != 1.0 or g != 0.0):
                        continue
                    o = simulate(S, 60000, 6.0, 0.3, d, ls, g, alb)
                    o.update(kind='scatter', g=g, ls=ls, d=d, alb=alb, P=6.0, s=0.3,
                             tau=d / ls, tau_tr=d * (1 - g) / ls)
                    rows.append(o)
                    print(f"g {g} ls {ls} d {d} alb {alb}: tau* {o['tau_tr']:.2f} open eta_beta {o['open_beta']:.3f} "
                          f"(x{o['open_beta'] / S['ref']:.2f}) crosstalk {o['crosstalk']:.2f} | bank Rw=.95 eta_beta "
                          f"{o['bank0.95_beta']:.3f} (x{o['bank0.95_beta'] / S['ref']:.2f}) Rw=.98 {o['bank0.98_beta']:.3f} "
                          f"[{time.time() - t0:.0f}s]", flush=True)
    # (3) pixel pitch / spacer study at a representative scatterer (g 0.6, ls 1, d 1)
    for P in (3.0, 4.0, 6.0, 8.0, 12.0, 20.0):
        for s in (0.1, 0.3, 1.0):
            o = simulate(S, 60000, P, s, 1.0, 1.0, 0.6, 0.999)
            o.update(kind='pixel', g=0.6, ls=1.0, d=1.0, alb=0.999, P=P, s=s)
            rows.append(o)
            print(f"P {P} s {s}: open eta_beta {o['open_beta']:.3f} crosstalk {o['crosstalk']:.2f} hits {o['mean_hits']:.1f} | "
                  f"bank .90 {o['bank0.9_beta']:.3f} .95 {o['bank0.95_beta']:.3f} .98 {o['bank0.98_beta']:.3f} "
                  f".995 {o['bank0.995_beta']:.3f} | TIR-air {o['tir_air_beta']:.3f} TIR-1.25 {o['tir_n125_beta']:.3f}  (ref {S['ref']:.3f}) [{time.time() - t0:.0f}s]", flush=True)
    json.dump(dict(ref=S['ref'], eta_sub=S['eta_sub'], validation=v, rows=rows), open('pixel_results.json', 'w'))
    print('done')
