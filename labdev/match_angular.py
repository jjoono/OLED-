"""Rank simulated angular profiles against measured ones.

Usage: python3 match_angular.py measured.csv [out.csv]
measured.csv columns: set (orange|green), structure (bare|freeform|hemisphere), device, trial, I0, I10, ..., I80
(I/I(0) on the measured side, 10-degree steps). Each measurement is compared with every simulated case of its set
and structure (and, for lens devices, every lens variant / patch size available); metric = rms deviation over
10-60 deg and 10-80 deg, plus the angular factor 2*int(I/I0) sin(t) dt (natural spline, I(90) = 0)."""
import sys, csv, json, glob, numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import trapezoid
A = np.arange(0, 90, 10); E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))


def rel(prof):
    p = np.array(prof); i0 = (p[:5] * dOm[:5]).sum() / dOm[:5].sum()
    return np.array([1.0] + [(p[a - 3:a + 3] * dOm[a - 3:a + 3]).sum() / dOm[a - 3:a + 3].sum() / i0 for a in A[1:]])


def pe(v):
    cs = CubicSpline(np.radians(np.r_[A, 90]), np.r_[v, 0.0], bc_type='natural'); t = np.linspace(0, np.pi / 2, 2001)
    return 2 * trapezoid(cs(t) * np.sin(t), t)


def sims():
    S = {'orange': {}, 'green': {}}
    R = json.load(open('layout_results.json'))
    for k, v in R.items():
        st, patch, rect = k.split('|')
        if rect == 'rectAl': S['orange'][(st, patch.replace('patch0mm', 'bare'))] = rel(v['profile'])
    for f in glob.glob('layout_variants_*.json'):
        for k, v in json.load(open(f)).items():
            S['orange'][('hemisphere' if k.startswith(('hemi', 'cap')) else 'freeform', 'patch3mm ' + k)] = rel(v['profile'])
    for f, tag in (('layout_green_results.json', 'ETL42'), ('layout_green_results_etl60.json', 'ETL60')):
        for k, v in json.load(open(f)).items():
            st, patch = k.split('|'); S['green'][(st, patch.replace('patch0mm', 'bare') + ' ' + tag)] = rel(v['profile'])
    return S


if __name__ == '__main__':
    S = sims(); rows = []
    for r in csv.DictReader(open(sys.argv[1])):
        m = np.array([float(r[f'I{a}']) for a in A]); best = []
        for (st, case), s in S[r['set']].items():
            if st != r['structure']: continue
            best.append((float(np.sqrt(np.mean((s[1:7] - m[1:7]) ** 2))), float(np.sqrt(np.mean((s[1:] - m[1:]) ** 2))), case, pe(s)))
        best.sort()
        for rk, (e60, e80, case, pes) in enumerate(best[:3]):
            rows.append([r['set'], r['structure'], r['device'], r['trial'], round(pe(m), 3), rk + 1, case, round(e60, 3), round(e80, 3), round(pes, 3)])
    w = csv.writer(open(sys.argv[2], 'w', newline='') if len(sys.argv) > 2 else sys.stdout)
    w.writerow(['set', 'structure', 'device', 'trial', 'PE_meas', 'rank', 'sim_case', 'rms_10-60', 'rms_10-80', 'PE_sim']); w.writerows(rows)
