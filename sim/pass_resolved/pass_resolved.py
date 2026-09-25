"""Pass-resolved escape probability p_k and round-trip loss A'_k of the matrix series
(Supplementary Note 1, eq. S7): eta_ext = sum_k p_k prod_{j<k} (1 - p_j)(1 - A'_j).
Fig. 2(d)-(f) generic stack at 550 nm (ITO 50 / organics 420 nm, n = 1.8), hemispherical MLA (n_MLA = n_sub)."""
import numpy as np, sys, csv, os
sys.path[:0] = ['/home/user/OLED-/sim/audit', '/home/user/OLED-/sim/fig3c']
import recompute_series as RS, cps2
H = os.path.dirname(os.path.abspath(__file__))
def passes(S, n, K=300):
    BT, BR, _ = RS.bsdf(n)
    R = np.real(cps2.stack_reflectance(S, RS.TH)); P = np.real(cps2.sub_angular(S, RS.TH))
    u = P / P.sum(); out = []
    for k in range(K):
        U = u.sum(); pk = (BT @ u) / U; w = BR @ u; Ak = 1 - (R @ w) / w.sum()
        out.append((k, U, pk, Ak, pk * U, np.sum(u[RS.TH > 60]) / U)); u = R * w
    tot = sum(o[4] for o in out)
    pc = np.sum(BT * RS.W) / np.sum(RS.W); Ac = 1 - np.sum(R * RS.W) / np.sum(RS.W)
    return out, tot, pc, Ac
rows = [['case', 'k', 'U_k (power at pass k)', 'p_k', "A'_k", 'eta_ext^(k)', 'fraction >60 deg', 'cumulative eta_ext', 'p (cos.sin)', "A' (cos.sin)", 'eq3']]
for name, refl in (('Ag', RS.M.AG), ('Al', RS.AL)):
    for n in (1.5, 1.8):
        out, tot, pc, Ac = passes(RS.stack(refl, n), n); c = 0
        for o in out:
            c += o[4]; rows.append([f'{name} n_sub {n}', *o, c, pc, Ac, pc / (pc + (1 - pc) * Ac)])
        print(f'{name} {n}: series {tot:.4f} eq3 {pc/(pc+(1-pc)*Ac):.4f} | p_k ' + ' '.join(f'{o[2]:.3f}' for o in out[:6]) + " | A'_k " + ' '.join(f'{o[3]:.4f}' for o in out[:6]) + f' | p {pc:.3f} A {Ac:.4f}')
with open(os.path.join(H, 'pass_resolved.csv'), 'w', newline='') as f: csv.writer(f).writerows(rows)
