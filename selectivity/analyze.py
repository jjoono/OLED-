"""Tables from the sweep: ceilings, best realistic DBRs, references.  Writes results/*.csv."""
import os, csv, json, numpy as np
import sel, common, bound

os.makedirs('results', exist_ok=True)
SW = common.Sweep()
MLA = dict(np.load('data/ref_mla.npz'))
CAV = np.load('data/ref_cavity.npz') if os.path.exists('data/ref_cavity.npz') else None
BETA_T = (10.0, 20.0, 30.0)


def wcsv(name, header, rows):
    with open(f'results/{name}.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)


def mla(stack, spec):
    return spec @ (SW.eta_sub[stack][:, None] * MLA[stack])


def cavity_best(spec, beta, cath, key='cone'):
    """Best strong-microcavity (thin Ag anode, planar) over the ETL/HTL scan."""
    if CAV is None:
        return None, None
    E = np.einsum('l,kl,klb->kb', spec, CAV['eta_sub'], CAV['cones'])
    sel_k = [i for i, k in enumerate(CAV['keys']) if k.startswith(cath + '|')]
    col = common.BI[5.0] if key == 'fwd' else common.BI[beta]
    i = max(sel_k, key=lambda i: E[i, col])
    return str(CAV['keys'][i]), E[i]


def fwd(e):
    return e[common.BI[5.0]] / np.sin(np.radians(5.0)) ** 2


# ---- 1. closed-form ceiling
rows = []
for n in (1.51, 1.8):
    for A in (0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.15, 0.2, 0.3):
        for b in (10, 20, 30, 45):
            rows.append([n, A, b, bound.eta_cone(A, b, n), bound.eta_open(A, n), bound.fwd_gain(A, b, n)])
wcsv('bound_closed_form', ['n_sub', 'A_prime', 'beta_deg', 'eta_beta_max', 'eta_open_no_filter', 'forward_luminance_gain'], rows)

# ---- 2. stack summary at each cone
rows = []
for st in sel.STACKS:
    for sp in ('G30', 'Ir(ppy)2acac', 'v-DABNA (HF)'):
        spec = common.SPECTRA[sp]
        m = mla(st, spec)
        planar = SW.ref(st, 'none', spec, ('fresnel',))
        for b in BETA_T:
            r = [st, sp, b, planar[common.BI[b]], m[common.BI[b]], fwd(planar), fwd(m)]
            for rnd in ('full', 'scat', 'none'):
                nof = SW.ref(st, rnd, spec, ('fresnel',))
                ide = SW.ref(st, rnd, spec, ('ideal', b))
                d, e = SW.best_edge(st, rnd, spec, b)
                r += [nof[common.BI[b]], ide[common.BI[b]], e[common.BI[b]], e[-1], fwd(nof), fwd(ide), fwd(e), '/'.join(map(str, d[1:]))]
            rows.append(r)
hdr = ['stack', 'spectrum', 'beta', 'planar_EQE_b', 'MLA_EQE_b', 'planar_fwd', 'MLA_fwd']
for rnd in ('full', 'scat', 'none'):
    hdr += [f'{rnd}_nofilter_EQE_b', f'{rnd}_ideal_EQE_b', f'{rnd}_DBR_EQE_b', f'{rnd}_DBR_EQE_total',
            f'{rnd}_nofilter_fwd', f'{rnd}_ideal_fwd', f'{rnd}_DBR_fwd', f'{rnd}_DBR_design']
wcsv('stack_summary', hdr, rows)

# ---- 3. bandwidth dependence (best DBR per material and pair count)
rows = []
for st in sel.STACKS:
    for rnd in ('full', 'scat'):
        for sp in common.SPECTRA:
            spec = common.SPECTRA[sp]
            for b in (10.0, 20.0):
                ide = SW.ref(st, rnd, spec, ('ideal', b))[common.BI[b]]
                nof = SW.ref(st, rnd, spec, ('fresnel',))[common.BI[b]]
                for mat in sel.DBR_MATS:
                    for N in (4, 6, 8, 12, 16):
                        d, e = SW.best_edge(st, rnd, spec, b, mat=mat, N=N)
                        rows.append([st, rnd, sp, common.FWHM[sp], b, mat, N, d[3], e[common.BI[b]], ide, nof, e[-1], fwd(e)])
wcsv('dbr_bandwidth', ['stack', 'randomizer', 'spectrum', 'fwhm_nm', 'beta', 'material', 'pairs', 'lam_edge_nm',
                       'EQE_beta', 'ideal_EQE_beta', 'nofilter_EQE_beta', 'EQE_total', 'EQE_fwd'], rows)

# ---- 4. strong-cavity reference
if CAV is not None:
    rows = []
    for sp in ('G30', 'Ir(ppy)2acac', 'v-DABNA (HF)'):
        spec = common.SPECTRA[sp]
        for cath in ('Al', 'Ag'):
            for b in BETA_T:
                k, e = cavity_best(spec, b, cath)
                rows.append([sp, cath, b, k, e[common.BI[b]], e[-1], fwd(e)])
            k, e = cavity_best(spec, 5.0, cath, key='fwd')
            rows.append([sp, cath, 'fwd', k, e[common.BI[5.0]], e[-1], fwd(e)])
    wcsv('strong_cavity', ['spectrum', 'cathode', 'objective', 'cathode|ETL|HTL', 'EQE_beta', 'EQE_total', 'EQE_fwd'], rows)

print(open('results/stack_summary.csv').read()[:4000])
