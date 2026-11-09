"""Shared post-processing: load the sweep, pick the best edge DBR for a spectrum and cone."""
import numpy as np, ast, sel

SPECTRA = {f'G{f:g}': sel.gauss(530.0, f) for f in (10, 20, 30, 45, 60, 80, 100)}
SPECTRA.update({k: sel.real_spec(k) for k in ('v-DABNA (HF)', 'blue QD', 'Ir(ppy)2acac', 'Ir(dmppy-ph)2tmd')})
FWHM = {f'G{f:g}': f for f in (10, 20, 30, 45, 60, 80, 100)}
FWHM.update({'v-DABNA (HF)': 23, 'blue QD': 22, 'Ir(ppy)2acac': 63, 'Ir(dmppy-ph)2tmd': 68})
BI = {b: i for i, b in enumerate(sel.BETAS)}


class Sweep:
    def __init__(self, path='data/sweep.npz'):
        z = np.load(path)
        self.designs = [tuple(float(x) if isinstance(x, np.floating) else x for x in eval(s, {'np': np})) for s in z['designs']]
        self.data = {tuple(k.split('|')): z[k] for k in z.files if '|' in k}
        DV = np.load('data/devices.npz')
        self.eta_sub = {k: DV[f'{k}_eta_sub'] for k in sel.STACKS}
        self.idx = {d: i for i, d in enumerate(self.designs)}

    def eqe(self, stack, rand, spec):
        """(n_design, n_beta) EQE into each cone for one spectrum (photon basis)."""
        c = self.data[(stack, rand)].astype(float)                 # (design, lam, beta)
        return np.einsum('l,dlb->db', spec * self.eta_sub[stack], c)

    def best_edge(self, stack, rand, spec, beta, mat=None, N=None, key='cone'):
        E = self.eqe(stack, rand, spec)
        col = BI[5.0] if key == 'fwd' else BI[beta]
        cand = [i for i, d in enumerate(self.designs) if d[0] == 'edge' and (mat is None or d[1] == mat) and (N is None or d[2] == N)]
        i = max(cand, key=lambda i: E[i, col])
        return self.designs[i], E[i]

    def ref(self, stack, rand, spec, design):
        return self.eqe(stack, rand, spec)[self.idx[design]]
