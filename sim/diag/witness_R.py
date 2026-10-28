"""Expected 30-deg s-pol reflectance (incidence from air on the glass side, n 1.52, front face incoherent) of the three
cathode-interface witnesses:  W1 glass | B3PyMPM:Cs2CO3 10 wt% d | Ag 100;  W2 ... | 50 wt% layer 5 nm | Ag;  W3 ... | Cs2CO3 2 nm | Ag.
Interlayers n 1.7 + i k; Ag: McPeak (ideal) or n + 0.1 (fit to the aux-reflector air-side data).  -> witness_R.csv"""
import numpy as np, csv
import refl_center_fit as C
D = C.D; NG = 1.52; A = np.radians(30.0)
ct = np.cos(np.arcsin(np.sin(A) / NG)); R1 = ((np.cos(A) - NG * ct) / (np.cos(A) + NG * ct)) ** 2
def R(lam, d_etl, inter, dag):
    j = int(round(lam - 400)); ns = [D.B3o[j]]; ds = [d_etl]
    if inter: ns.append(1.7 + 1j * inter[1]); ds.append(inter[0])
    ns += [D.AG[j] + dag, 1.0]; ds += [100.0, 0]
    r = C.tmm_s(ns, ds, lam, NG, 2 * np.pi / lam * np.sin(A)); R2 = abs(r) ** 2
    return R1 + (1 - R1) ** 2 * R2 / (1 - R1 * R2)
CASES = [('W1 doped/Ag', None), ('W2 +50wt% 5nm k0', (5, 0.0)), ('W2 +50wt% 5nm k0.3', (5, 0.3)),
         ('W3 +Cs2CO3 2nm k0', (2, 0.0)), ('W3 +Cs2CO3 2nm k0.3', (2, 0.3)), ('W3 +Cs2CO3 2nm k1', (2, 1.0)), ('W3 +Cs2CO3 2nm k2', (2, 2.0))]
L = np.arange(450, 701, 1); P = (450, 500, 550, 600, 650, 700)
rows = []
for d in (60.0,):
    for dag, an in ((0.0, 'ideal Ag'), (0.1, 'Ag n+0.1')):
        print('--- doped layer %d nm, %s ---' % (d, an))
        for name, inter in CASES:
            r = np.array([R(l, d, inter, dag) for l in L])
            row = dict(d_etl=d, Ag=an, case=name, **{'R%d' % p: round(100 * r[p - 450], 1) for p in P}, mean=round(100 * r.mean(), 1))
            rows.append(row); print('%-22s ' % name + '  '.join('%d:%.1f' % (p, 100 * r[p - 450]) for p in P) + '  mean %.1f' % (100 * r.mean()))
with open('witness_R.csv', 'w', newline='') as fo:
    w = csv.DictWriter(fo, list(rows[0])); w.writeheader(); w.writerows(rows)
