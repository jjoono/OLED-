"""Builds supplementary Note S12 (technical-review robustness checks) from data/review_checks.json and data/review_traces.json."""
import json
C = json.load(open('data/review_checks.json')); T = json.load(open('data/review_traces.json'))
L = []
A = L.append
A('## Note S12 | Device inputs and robustness of the spectral metrics\n')
A('**Device inputs.** Table S14 lists the three stacks; Fig. S11 shows the first-pass radiance launched into the substrate and the angle-resolved reflectance at 550 nm. The emitters form a single dipole sheet at the centre of the 25 nm birefringent EML (ordinary and extraordinary indices from measured data), with a horizontal dipole ratio of 0.865. Dipole emission and reflectance are computed for s and p polarization; $R(\\theta)$ is their average and $\\mathbf{v}^{(0)}$ their sum.\n')
A('**Table S14 | Device stacks (from the cathode to the substrate).**\n')
A('| Device | Cathode | ETL | EML | HTL-side layers | Anode | Substrate | $A\'$ at 550 nm |')
A('|---|---|---|---|---|---|---|---|')
A('| Al | Al 100 nm | B3PyMPM 61.3 nm | TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd 25 nm | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.167 |')
A('| Ag | Ag 100 nm | B3PyMPM 61.3 nm | as Al | TCTA 10 nm / TAPC 97.5 nm | ITO 150 nm | glass, $n$ = 1.51 | 0.044 |')
A('| Low loss | Ag 100 nm | B3PyMPM 200 nm | as Al | TAPC 180 nm | ITO 50 nm | glass, $n$ = 1.51 | 0.026 |\n')
A('**Fig. S11 | Device inputs at 550 nm.** (a) Radiance launched into the substrate on the first pass, $v^{(0)}_k/w_k$, relative to its mean over 0–20°. (b) Angle-resolved reflectance $R(\\theta)$ of the stack seen from the substrate, with the uniform-radiance round-trip loss $A\'$. The dotted line marks the critical angle of the glass/air interface. *(`figS11_device_inputs.png`)*\n')
A('**Discretization, cutoff, ray count, and domain size.** Table S15 reports $N_{\\mathrm{mix}}$, $\\xi$, $\\eta_{\\mathrm{sta}}$ and the reciprocal reduced-model error (Al, 550 nm) for representative structures under changes of the polar bin width, of the cutoff angle of the spectral metrics, of the number of rays per incidence bin, and of the size of the periodic simulation domain. Bin widths above 1° are obtained by étendue-weighted merging; the scattering layer is also evaluated on its native 0.5° grid. $\\eta_{\\mathrm{sta}}$ is insensitive to all of these choices. $N_{\\mathrm{mix}}$ and $\\xi$ change with the cutoff because the near-grazing states, which carry little étendue, can retain angle; the ranking of the structures is preserved (rank correlation with the 85° values: $N_{\\mathrm{mix}}$ 0.959, 0.991, 0.962 and $\\xi$ 0.902, 0.979, 0.908 for 80°, 88°, 90°).\n')
A('**Table S15 | Convergence of the metrics (Al, 550 nm).**\n')
A('| Structure | Variant | $N_{\\mathrm{mix}}$ | $\\xi$ | $\\eta_{\\mathrm{sta}}$ | reduced − full (pp) |')
A('|---|---|---|---|---|---|')
lab = {'1deg_cut80': '1° bins, cutoff 80°', '1deg_cut85': '1° bins, cutoff 85° (reference)', '1deg_cut88': '1° bins, cutoff 88°', '1deg_cut90': '1° bins, no cutoff', '2deg_cut85': '2° bins, cutoff 85°', '3deg_cut85': '3° bins, cutoff 85°', '0.5deg_cut85': '0.5° bins (native), cutoff 85°'}
nice = {'planar': 'planar', 'MLA hemi_total': 'hemispherical MLA', 'scatter g0.8_t1.0': 'scattering layer ($g$ 0.8, $\\tau^*$ 1)', 'A_shape1.0': 'random lenslet shapes, $s$ = 1', 'W_slope0.2': 'wavy, rms slope 0.2', 'W_slope1.0': 'wavy, rms slope 1.0'}
fmtN = lambda v: '∞' if v > 1e6 else f'{v:.2f}'
for nm, d in C['convergence'].items():
    for k in ('0.5deg_cut85', '1deg_cut85', '2deg_cut85', '3deg_cut85', '1deg_cut80', '1deg_cut88', '1deg_cut90'):
        if k in d:
            v = d[k]; A(f"| {nice[nm]} | {lab[k]} | {fmtN(v['Nmix'])} | {v['xi']:.3f} | {v['eta']:.4f} | {v['err_pp']:+.2f} |")
rows = [('A_shape0.0_N5000', 'hemispherical array ($s$ = 0)', '5,000 rays per bin'), ('A_shape0.0_N10000', 'hemispherical array ($s$ = 0)', '10,000 rays per bin'), ('A_shape0.0_library', 'hemispherical array ($s$ = 0)', '20,000 rays per bin (library)'),
        ('A_shape0.0_N40000', 'hemispherical array ($s$ = 0)', '40,000 rays per bin'), ('W_slope0.2_N2500', 'wavy, rms slope 0.2', '2,500 rays per bin'), ('W_slope0.2_N5000', 'wavy, rms slope 0.2', '5,000 rays per bin'),
        ('W_slope0.2_library', 'wavy, rms slope 0.2', '10,000 rays per bin (library)'), ('W_slope0.2_N20000', 'wavy, rms slope 0.2', '20,000 rays per bin'),
        ('A_shape1.0_library', 'random lenslet shapes, $s$ = 1', '6 × 6 supercell (library)'), ('A_shape1.0_K10', 'random lenslet shapes, $s$ = 1', '10 × 10 supercell, new realization'),
        ('C_jitter0.3_library', 'position jitter 0.3', '6 × 6 supercell (library)'), ('C_jitter0.3_K10', 'position jitter 0.3', '10 × 10 supercell, new realization'),
        ('W_slope0.2_library', 'wavy, rms slope 0.2', '12 × 12 box (library)'), ('W_slope0.2_L24', 'wavy, rms slope 0.2', '24 × 24 box, new realization'),
        ('W_slope1.0_library', 'wavy, rms slope 1.0', '12 × 12 box (library)'), ('W_slope1.0_L24', 'wavy, rms slope 1.0', '24 × 24 box, new realization')]
for k, s_, v_ in rows:
    if k in T:
        v = T[k]; A(f"| {s_} | {v_} | {v['Nmix']:.2f} | {v['xi']:.3f} | {v['eta']:.4f} | {v['err_pp']:+.2f} |")
A('')
S = C['symmetrized']; B = C['bounds']
A(f"**Detailed balance and bounds.** Symmetrizing $F_{{ji}} = w_i B_R(i\\to j)$ to $(F+F^{{\\mathsf T}})/2$ and rescaling each column to its traced return fraction, which keeps $p(\\theta)$, changes $N_{{\\mathrm{{mix}}}}$ by at most {S['max_dNmix']:.2f} ({100*S['max_rel_dNmix']:.1f}%), $\\xi$ by at most {S['max_dxi']:.3f}, $\\eta_{{\\mathrm{{sta}}}}$ by at most {S['max_deta_pp']:.2f} pp, and the reduced-model error by at most {S['max_derr_pp']:.2f} pp over all 37 non-planar structures and three devices. The bounds $\\min_\\theta p(\\theta) \\le p_\\infty \\le \\tilde p$ hold for every structure, with smallest margins of 0.0014 (upper) and 0.106 (lower, thick scattering layer); the planar interface meets the lower bound with equality ($p_\\infty$ = 0). With $p$ in place of $\\tilde p$ the upper bound would fail for the thick scattering layer ($\\tau^*$ = 4) by 0.005, because return beyond the cutoff then counts as neither escape nor return. The device-aware variant $1-\\rho(\\mathrm{{diag}}\\,\\mathbf{{R}}\\,B_R)$ (Al, 550 nm) ranks the structures like $p_\\infty$ (Spearman {C['device_aware']['spearman']:.3f}). The full-series solve is well conditioned: the largest spectral radius of $\\mathrm{{diag}}\\,\\mathbf{{R}}\\,B_R$ is {C['solver']['rho_max']:.4f}, the largest residual {C['solver']['residual_max']:.0e}, and a pass-by-pass sum agrees with the solve to {C['solver']['sum_vs_solve_max']:.0e}.\n")
A('**Contour templates.** Fig. 3a scales the angular shape of $p(\\theta)$ of the hemispherical MLA. Table S16 repeats the calculation with the shapes of a scattering layer and of a wavy surface (Al, 550 nm).\n')
A('**Table S16 | Sensitivity of $\\eta_{\\mathrm{sta}}$ to mixing and to escape for three $p(\\theta)$ templates.** Span: change of $\\eta_{\\mathrm{sta}}$ when $N_{\\mathrm{mix}}$ varies from 1 to 16 at fixed $p$. Slope: change of $\\eta_{\\mathrm{sta}}$ per 0.01 in $p$ at $N_{\\mathrm{mix}}$ = 4.\n')
A('| Template | span at $p$ = 0.30 (pp) | at 0.38 | at 0.42 | slope per 0.01 in $p$ (pp) |'); A('|---|---|---|---|---|')
tn = {'MLA hemi_total': 'hemispherical MLA', 'scatter g0.8_t1.0': 'scattering layer ($g$ 0.8, $\\tau^*$ 1)', 'W_slope1.0': 'wavy surface, rms slope 1.0'}
for nm, v in C['templates'].items():
    A(f"| {tn[nm]} | {v['rows'][0]['span_pp']:.2f} | {v['rows'][1]['span_pp']:.2f} | {v['rows'][2]['span_pp']:.2f} | {v['deta_dp_per_0p01_pp']:.2f} |")
A('')
A('**Counterfactual decomposition.** Along the lenslet-shape sweep, the change of $\\eta_{\\mathrm{sta}}$ relative to the hemispherical array ($s$ = 0) is split in the reciprocal reduced model into a mixing part (mixing $m$ of the randomized array, $p(\\theta)$ of the hemisphere), an escape-shape part ($p(\\theta)$ of the randomized array, $m$ of the hemisphere), and their interaction (Table S17).\n')
A('**Table S17 | Decomposition of the change in $\\eta_{\\mathrm{sta}}$ along the shape sweep (550 nm, pp).**\n')
A('| $s$ | Device | mixing only | $p(\\theta)$ only | interaction | total (reduced) | total (full series) |'); A('|---|---|---|---|---|---|---|')
for r in C['decomposition']:
    A(f"| {r['s']:.1f} | {'low loss' if r['dev']=='LL15' else r['dev']} | {r['mixing']:+.2f} | {r['shape']:+.2f} | {r['interaction']:+.2f} | {r['total']:+.2f} | {r['full']:+.2f} |")
A('')
A('**Composition window.** Table S18 gives the range of the selectivity $S_{40\\text{–}60} = P_{\\mathrm{air}}(40\\text{–}60°)/P_{\\mathrm{air}}(0\\text{–}90°)$ over all 38 structures for each device and wavelength.\n')
A('**Table S18 | Range of $S_{40\\text{–}60}$ over all structures.**\n'); A('| Device, wavelength | min | max |'); A('|---|---|---|')
dn = {'Al|500.0': 'Al, 500 nm', 'Al|550.0': 'Al, 550 nm', 'Al|600.0': 'Al, 600 nm', 'Ag|550.0': 'Ag, 550 nm', 'LL15|550.0': 'low loss, 550 nm'}
for k, (a, b) in C['window'].items(): A(f'| {dn[k]} | {a:.3f} | {b:.3f} |')
A('\nIn the low-loss device the lower end (0.274–0.293) belongs to the hemispherical MLAs; all other structures stay within 0.31–0.37.\n')
open('../note_s12.md', 'w').write('\n'.join(L))
print('note written', len(L))
