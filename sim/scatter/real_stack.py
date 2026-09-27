"""Real-stack comparison at 550 nm: Ag 100 / ETL d_ETL / TCTA:B3PyMPM 25 (centre) / TAPC d_HTL / ITO 50 (k 0.0032) / n_sub 1.8.
ETL: B3PyMPM (n_o 1.82, n_e 1.61) or 3TPYMB (isotropic 1.61); Theta 0.67 or 0.9.
Extraction structures (all n 1.8): hemisphere MLA; MLA on a weak scattering layer; scattering layer alone."""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/scatter']
import cps2, materials as M
from mla_plus_scatter import slab, BTm, BRm, TH, g0
m = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)['material']
TY = complex(np.asarray(m.TPYMB_3)[150])
EML = (1.83273, 1.67122); TAPC = (1.69075, 1.66375); B3 = M.ETL['B3PyMPM']
def combo(S, g):
    T, R = slab(S, g, N=8000, seed=int(S*97 + g*13)); X = np.linalg.solve(np.eye(90) - R @ BRm, T)
    return BTm @ X, R + T @ BRm @ X
scat = np.load('/home/user/OLED-/sim/scatter/scatter_bsdf_n1.8.npz')
STRUCT = {'MLA': (BTm, BRm), 'MLA+scat(S1,g0.9)': combo(1, 0.9), 'MLA+scat(S0.3,g0.5)': combo(0.3, 0.5),
          'scat only (S10,g0.9)': (scat['S10_g0.9'][:90].sum(0), scat['S10_g0.9'][90:][::-1])}
def run(a):
    etl, th, de, dh = a; ne = (B3 if etl == 'B3PyMPM' else (TY, TY))
    S = cps2.Stack(550.0, EML, 25.0, 12.5, above=[(ne[0], ne[1], de), (M.AG, M.AG, 100.0)],
                   below=[(TAPC[0], TAPC[1], dh), (M.ITO, M.ITO, 50.0)], n_sub=1.8, h=th)
    r = cps2.solve_pol(S, npts=8000); es = float(np.real(r['air'] + r['sub']))
    R = np.real(cps2.stack_reflectance(S, TH)); P = np.real(cps2.sub_angular(S, TH)); P /= P.sum(); A = 1 - R @ g0
    out = [etl, th, de, dh, es, A, float(P[TH > 70].sum())]
    for name, (BT, BR) in STRUCT.items():
        v = P.copy(); e = 0
        for k in range(400): e += BT @ v; v = R*(BR @ v)
        pl = BT @ g0; e3 = pl/(pl + (1 - pl)*A); out += [e, e3, es*e]
    return out
if __name__ == '__main__':
    jobs = [(etl, th, de, dh) for etl in ('B3PyMPM', '3TPYMB') for th in (0.67, 0.9) for de in (200.0, 300.0, 400.0) for dh in (200.0, 300.0)]
    with Pool(4) as p: res = p.map(run, jobs)
    hdr = ['ETL', 'Theta', 'd_ETL', 'd_HTL', 'eta_sub', "A'", 'Psub>70']
    for n in STRUCT: hdr += [f'{n} eta_ext', f'{n} eq3', f'{n} EQE']
    with open('real_stack.csv', 'w', newline='') as f: csv.writer(f).writerows([hdr] + res)
    print('ETL      Th   dE  dH | eta_sub  A%  >70 | ' + ' | '.join(f'{n}: ext/EQE' for n in STRUCT))
    for r in res:
        print(f"{r[0]:8s} {r[1]:.2f} {r[2]:.0f} {r[3]:.0f} | {r[4]:.3f} {100*r[5]:.1f} {r[6]:.2f} | " + ' | '.join(f'{r[7+3*i]:.3f}/{r[9+3*i]:.3f}' for i in range(4)) + f"  (eq3 MLA {r[8]:.3f})")
