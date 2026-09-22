#!/usr/bin/env python3
"""Recompute every panel that was evaluated with the closed form of eq. (3) using the matrix
series of eq. (1) and the author's LightTools hemisphere BSDF (lt_hemisphere_bsdf.mat,
slices n_MLA = 1.30:0.05:2.00, n_MLA = n_sub), so that Figs 2(c), 2(d)-(f), 3(a), 5(c) and
Supplementary Table 1 sit on the same footing as Figs 4 and 5(a),(b),(d),(e).

Stack (550 nm, PLQY = 1, Theta = 2/3 unless swept):
    air | reflector 100 nm | ETL n=1.8 | EML n=1.8 20 nm, dipole at the centre | HTL n=1.8 | TCO | substrate n_sub
"""
import numpy as np, scipy.io as sio, csv, sys, os, time
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, '/home/user/OLED-/sim/mla')
import cps2, materials as M, series

HERE = os.path.dirname(os.path.abspath(__file__))
AL = 0.958336812040867 + 6.68678039868499j
TH = np.arange(90) + 0.5
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
B = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA']
NS = np.round(np.arange(1.30, 2.001, 0.05), 2)


def bsdf(n):
    k = int(np.argmin(abs(NS - n)))
    b = B[:, :, k]
    return b[:90, :].sum(0), b[90:180, :][::-1, :], NS[k]


def evaluate(S, n_mla, npts=12000):
    """eta_sub, A', p(cos.sin), p(P_sub), eta_ext(series), eta_ext(closed, p cos.sin), EQE(series), EQE(closed)."""
    BT, BR, _ = bsdf(n_mla)
    r = cps2.solve_pol(S, npts=npts)
    es = r['air'] + r['sub']
    R = cps2.stack_reflectance(S, TH); P = cps2.sub_angular(S, TH)
    A = 1 - np.sum(R * W) / np.sum(W)
    pc = np.sum(BT * W) / np.sum(W); pp = BT @ (P / P.sum())
    e, terms = series.eta_ext(BT, BR, R, P, n_term=100)
    ecf = pc / (pc + (1 - pc) * A)
    return dict(eta_sub=es, spp=r['spp'], wg=r['wg'], abs=r['abs'], Aprime=A, p_cos_sin=pc, p_Psub=pp,
                eta_ext_series=e, eta_ext_closed=ecf, EQE_series=es * e, EQE_closed=es * ecf,
                tail=terms[-1] / e)


def stack(refl, n_sub, d_etl=200.0, d_htl=200.0, tco=M.ITO, d_tco=50.0, h=2/3, eml=(1.8, 1.8), d_eml=20.0):
    return cps2.Stack(550.0, eml, d_eml, d_eml / 2, above=[(1.8, 1.8, d_etl), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, d_htl), (tco, tco, d_tco)], n_sub=n_sub, h=h)


def write(name, header, rows):
    with open(os.path.join(HERE, name), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header)
        for r in rows: w.writerow(r)
    print('  wrote', name, len(rows), 'rows', file=sys.stderr, flush=True)


KEYS = ['eta_sub', 'spp', 'wg', 'abs', 'Aprime', 'p_cos_sin', 'p_Psub', 'eta_ext_series', 'eta_ext_closed', 'EQE_series', 'EQE_closed', 'tail']
t0 = time.time()

# (1) Fig. 2(c): ITO thickness at n_sub = 1.5, Al and Ag  (the Fig. 2 stack has 150 nm ITO as the anode)
rows = []
for name, refl in (('Al', AL), ('Ag', M.AG)):
    for d in np.arange(30.0, 200.1, 10.0):
        v = evaluate(stack(refl, 1.5, d_tco=d), 1.5)
        rows.append([name, '%.0f' % d] + ['%.6f' % v[k] for k in KEYS])
        print('  fig2c %s ITO %3.0f  eta_ext series %.4f closed %.4f' % (name, d, v['eta_ext_series'], v['eta_ext_closed']), file=sys.stderr, flush=True)
write('fig2c_series.csv', ['reflector', 'd_ITO_nm'] + KEYS, rows)

# (2) Fig. 2(d)-(f): substrate / MLA index, ITO 50 nm (as in the panel as drawn)
rows = []
for name, refl in (('Al', AL), ('Ag', M.AG)):
    for n in NS:
        v = evaluate(stack(refl, float(n), d_tco=50.0), float(n))
        rows.append([name, '%.2f' % n] + ['%.6f' % v[k] for k in KEYS])
        print('  fig2def %s n %.2f  p %.4f  eta_ext series %.4f closed %.4f  EQE series %.4f closed %.4f' % (name, n, v['p_cos_sin'], v['eta_ext_series'], v['eta_ext_closed'], v['EQE_series'], v['EQE_closed']), file=sys.stderr, flush=True)
write('fig2def_series.csv', ['reflector', 'n_sub'] + KEYS, rows)

# (3) Fig. 3(a): TCO extinction coefficient (ITO 50 nm, n = 1.8636) and thin-Ag real index (10 nm, k = 3.819), n_sub = 1.8, Ag reflector
rows = []
for k_ito in np.arange(0.0, 0.0801, 0.005):
    v = evaluate(stack(M.AG, 1.8, tco=complex(M.ITO.real, k_ito)), 1.8)
    rows.append(['ITO_k', '%.4f' % k_ito] + ['%.6f' % v[x] for x in KEYS])
    print('  fig3a k_ITO %.4f  EQE series %.4f closed %.4f' % (k_ito, v['EQE_series'], v['EQE_closed']), file=sys.stderr, flush=True)
for n_ag in np.arange(0.0, 0.501, 0.025):
    tco = complex(n_ag, 3.818978372838121)
    v = evaluate(stack(M.AG, 1.8, tco=tco, d_tco=10.0), 1.8)
    rows.append(['Ag_n', '%.4f' % n_ag] + ['%.6f' % v[x] for x in KEYS])
    print('  fig3a n_Ag %.3f  EQE series %.4f closed %.4f' % (n_ag, v['EQE_series'], v['EQE_closed']), file=sys.stderr, flush=True)
write('fig3a_series.csv', ['sweep', 'value'] + KEYS, rows)

# (4) Fig. 5(c): dipole orientation.  reference Al, n 1.5, 80/230 nm (cavity optimum), ITO 50; proposed Ag, n 1.8, 200/200
rows = []
for dev, refl, n, de, dh in (('reference', AL, 1.5, 80.0, 230.0), ('proposed', M.AG, 1.8, 200.0, 200.0)):
    for h in np.arange(0.5, 1.001, 0.05):
        v = evaluate(stack(refl, n, d_etl=de, d_htl=dh, h=float(h)), n)
        rows.append([dev, '%.2f' % h] + ['%.6f' % v[x] for x in KEYS])
        print('  fig5c %s Theta %.2f  EQE series %.4f closed %.4f' % (dev, h, v['EQE_series'], v['EQE_closed']), file=sys.stderr, flush=True)
write('fig5c_series.csv', ['device', 'Theta'] + KEYS, rows)

# (5) Supplementary Table 1: the real material stack, n_sub 1.8, hemisphere MLA (slice 11)
B3 = M.ETL['B3PyMPM']; TAPC = (1.69075, 1.66375); EML = (1.83273, 1.67122)
rows = []
for d_etl in (200.0, 300.0):
    for k_ito in (0.0032, 0.002):
        for h in (2/3, 0.8, 0.9):
            ito = complex(M.ITO.real, k_ito)
            S = cps2.Stack(550.0, EML, 25.0, 12.5, above=[(B3[0], B3[1], d_etl), (M.AG, M.AG, 100.0)],
                           below=[(TAPC[0], TAPC[1], 180.0), (ito, ito, 50.0)], n_sub=1.8, h=h)
            v = evaluate(S, 1.8)
            rows.append(['%.0f' % d_etl, '%.4f' % k_ito, '%.3f' % h] + ['%.6f' % v[x] for x in KEYS])
            print('  table1 ETL %.0f k %.4f Theta %.2f  EQE series %.4f closed %.4f' % (d_etl, k_ito, h, v['EQE_series'], v['EQE_closed']), file=sys.stderr, flush=True)
write('table1_series.csv', ['d_ETL_nm', 'k_ITO', 'Theta'] + KEYS, rows)
print('done in %.0f s' % (time.time() - t0), file=sys.stderr, flush=True)
print('done')
