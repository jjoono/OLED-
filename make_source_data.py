"""Merge the per-figure source-data workbooks into one file, one sheet per panel.

Output: Figure_source_data.xlsx (README sheet first). Run after selectivity/oe_figs.py,
make_figures.py and make_optimizer_convergence.py have written the per-figure files.
"""
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

OUT = 'Figure_source_data.xlsx'
HDR = Font(bold=True); FILL = PatternFill('solid', fgColor='E8EEF6')

def rows_of(path, sheet):
    return [list(r) for r in load_workbook(path, read_only=True, data_only=True)[sheet].iter_rows(values_only=True)]

def pick(rows, cols):
    h = rows[0]; idx = [h.index(c) for c in cols]
    return [[r[i] for i in idx] for r in rows]

wb = Workbook(); readme = wb.active; readme.title = 'README'
index = []  # (sheet, figure, panel, description, source file)

def add(name, fig, panel, desc, rows, src, note=None):
    assert len(name) <= 31, name
    ws = wb.create_sheet(name); r0 = 1
    if note:
        ws.cell(1, 1, note).font = Font(italic=True); r0 = 3
    for i, row in enumerate(rows):
        for j, v in enumerate(row):
            c = ws.cell(r0 + i, j + 1, v)
            if i == 0:
                c.font = HDR; c.fill = FILL
    ws.freeze_panes = ws.cell(r0 + 1, 1)
    for j in range(min(len(rows[0]) if rows else 0, 40)):
        ws.column_dimensions[get_column_letter(j + 1)].width = 14 if j else 22
    index.append((name, fig, panel, desc, src))

D = 'figures_oe/'
# ---------------- Fig. 1
f = D + 'fig1_mixing_concept.xlsx'
add('Fig1a', 'Fig. 1', 'a', 'Schematic; no plotted data.', [['note'], ['Panel (a) is a schematic.']], f)
for sh, nm, lab in (('planar M', 'Fig1b_planar_M', 'planar interface'), ('hemispherical MLA M', 'Fig1c_hemiMLA_M', 'hemispherical MLA'),
                    ('scattering layer M', 'Fig1d_scatter_M', 'scattering layer (g 0.8, tau* 1)')):
    add(nm, 'Fig. 1', nm[4], f'Return matrix M of the {lab}; row = returned 1-deg bin, column = incident 1-deg bin (substrate angle). Plotted as log10 with floor 1e-3.', rows_of(f, sh), f)
# ---------------- Fig. 2
f = D + 'fig2_complete_mixing.xlsx'
add('Fig2a_p_mix', 'Fig. 2', 'a', 'Uniform-radiance escape p, its etendue-weighted s.d. sigma_p, and completely mixed escape p_mix = p - sigma_p^2/(1-p), all structures.', rows_of(f, 'a_p_mix'), f)
add('Fig2b_gain', 'Fig. 2', 'b', 'Gain over the actual eta_sta predicted by complete reciprocal mixing and by the fully randomized model, vs N_mix (Al, 550 nm).', rows_of(f, 'b_gain'), f)
add('Fig2c_placement', 'Fig. 2', 'c', 'Complete reciprocal mixing at the escape surface versus at the device for all structures (Al, 550 nm); last column is the plotted difference.', rows_of(f, 'c_placement'), f)
add('Fig2d_slab', 'Fig. 2', 'd', 'Slab of index n with an isotropic internal source, a planar exit face, and an angle-independent back reflectance R.', rows_of(f, 'd_slab'), f)
# ---------------- Fig. 3 (split one table into panels)
f = D + 'fig3_reduced_model.xlsx'; C = rows_of(f, 'cases'); key = ['structure', 'class', 'device', 'lambda_nm', 'eta_full']
for nm, fig, col, desc in (('FigS2a_fully_randomized', 'Fig. S2', 'eta_fully_randomized', 'fully randomized model'), ('FigS2b_angle_preserving', 'Fig. S2', 'eta_m0', 'angle-preserving limit (m = 0)'),
                      ('FigS2c_reduced_isotropic', 'Fig. S2', 'eta_reduced_w', 'reduced model, m = 1 - lambda2, isotropic q ~ w'),
                      ('Fig3a_reduced_reciprocal', 'Fig. 3', 'eta_reduced_recip', 'reduced model, m = 1 - lambda2, reciprocal q ~ w(1-p)')):
    add(nm, fig, nm.split('_')[0][-1], f'x = eta_full (full series), y = {col} ({desc}); 342 cases.', pick(C, key + [col]), f)
h = C[0]; i_f, i_r = h.index('eta_full'), h.index('eta_reduced_recip')
E = [['structure', 'class', 'device', 'lambda_nm', 'error_pp (reciprocal model - full)']] + [r[:4] + [(r[i_r] - r[i_f]) * 100] for r in C[1:]]
groups = {}
for r in E[1:]:
    groups.setdefault((r[2], r[3]), []).append(abs(r[4]))
E2 = [['device', 'lambda_nm', 'mean_abs_error_pp (bar)', 'max_abs_error_pp', 'n']] + [[d, l, sum(v) / len(v), max(v), len(v)] for (d, l), v in groups.items()]
add('Fig3b_error_points', 'Fig. 3', 'b', 'Signed error of the reciprocal reduced model (pp), per case (points).', E, f)
add('Fig3b_error_bars', 'Fig. 3', 'b', 'Mean absolute error per device and wavelength (bars).', E2, f)
add('Fig3_all_columns', 'Fig. 3', 'a-b, S2', 'Full case table incl. p_Lamb, lambda2 (theta<85 deg) and p_inf.', C, f)
# ---------------- Fig. 4
f = D + 'fig4_structure_map.xlsx'; S = rows_of(f, 'structures')
add('Fig4a_structures', 'Fig. 4', 'a', 'Points: p (uniform radiance) and N_mix = 1/(1-lambda2) (Al, 550 nm). N_mix plotted capped at 25; planar has N_mix = inf (stored as 1e9).', S, f)
add('Fig4a_contours', 'Fig. 4', 'a', 'Contours: eta_sta of the reciprocal reduced model on a (p, N_mix) grid; first column N_mix, first row p.', rows_of(f, 'contours'), f)
h = S[0]; iN = h.index('N_mix'); B = [['series', 'parameter', 'parameter / max', 'N_mix', 'structure']]
hemi0 = next(r for r in S if r[0] == 'A_shape0.0')[iN]
for pre, lab, k in (('A_shape', 'lenslet shape spread s', 'shape'), ('C_jitter', 'position jitter / pitch', 'jitter'),
                    ('W_slope', 'wavy surface, rms slope', 'slope'), ('B_aspect', 'height spread', 'aspect')):
    pts = sorted((float(r[0].split(k)[1].split('_')[0]), r[iN], r[0]) for r in S[1:] if r[0].startswith(pre))
    if pre in ('C_jitter', 'B_aspect'):
        pts = [(0.0, hemi0, 'A_shape0.0')] + pts
    xm = max(p[0] for p in pts)
    B += [[lab, x, x / xm, n, s] for x, n, s in pts]
ref = {r[0]: r[iN] for r in S[1:]}
B += [['reference line: hemispherical MLA', None, None, ref['MLA hemi_total'], 'MLA hemi_total'],
      ['reference line: scattering layer', None, None, ref['scatter g0.8_t1.0'], 'scatter g0.8_t1.0']]
add('Fig4b_randomness', 'Fig. 4', 'b', 'N_mix against each randomness parameter normalized by its maximum (Al, 550 nm); dotted reference lines at the end.', B, f)
add('Fig4c_xi', 'Fig. 4', 'c', 'Late-pass escape ratio xi = p_inf / p~ against N_mix (Al, 550 nm).', pick(S, ['structure', 'class', 'N_mix', 'p_inf', 'xi']), f)
# ---------------- Fig. 5
f = D + 'fig5_bound_window.xlsx'
add('FigS4a_first_pass', 'Fig. S4', 'a', 'First-pass power into each polar band / substrate power (Al, 550 nm): bounds and four structures.', rows_of(f, 'a_first_pass'), f)
add('FigS4b_composition', 'Fig. S4', 'b', 'S40-60 of the multi-pass output vs eta_sta, all structures (Al, 550 nm); shaded band 0.33-0.36, Lambertian 0.337.', rows_of(f, 'b_composition'), f)
# ---------------- Fig. 6
f = D + 'fig6_hemisphere_benchmark.xlsx'
add('Fig5a_profiles', 'Fig. 5', 'a', 'Band-optimized freeform profiles and the hemispherical reference (units of lens radius).', rows_of(f, 'a_profiles'), f)
add('Fig5b_best_EQE', 'Fig. 5', 'b', 'Best EQE of each lens class per objective and G_j.', rows_of(f, 'b_best_EQE'), f)
add('Fig5c_searches', 'Fig. 5', 'c', 'Freeform/hemisphere ratio of each independent search.', rows_of(f, 'c_independent_searches'), f)
# ---------------- Supplement
f = D + 'figS11_device_inputs.xlsx'
add('FigS1_device_inputs', 'Fig. S1', 'a-b', 'First-pass substrate radiance (a) and angle-resolved reflectance R (b) of the three devices at 550 nm.', rows_of(f, 'devices_550nm'), f)
f = D + 'figS2_composition_window.xlsx'
add('FigS3_composition', 'Fig. S3', 'all', 'S40-60 vs eta_sta for all structures, devices and wavelengths.', rows_of(f, 'points'), f)
f = 'fig1_platform.xlsx'
add('FigS5b_profiles', 'Fig. S5', 'b', 'Hemispherical and freeform (40-60 deg) lens profiles.', rows_of(f, 'b_lens_profiles'), f)
add('FigS5b_control_pts', 'Fig. S5', 'b', 'Free control points of the freeform profile.', rows_of(f, 'b_control_points'), f)
add('FigS5c_bands', 'Fig. S5', 'c', 'Polar-band definition.', rows_of(f, 'c_band_definition'), f)
f = 'figS1_restart_control.xlsx'
add('FigS7a_restart_arms', 'Fig. S7', 'a', 'Re-measured hemisphere and restarted optimum for every arm.', rows_of(f, 'a_all_arms'), f)
add('FigS7b_confirm_reps', 'Fig. S7', 'b', 'Five-repeat confirmation of the 20-40 deg arm (repeats).', rows_of(f, 'b_confirmation_reps'), f)
add('FigS7b_confirm_stats', 'Fig. S7', 'b', 'Five-repeat confirmation (statistics).', rows_of(f, 'b_confirmation_stats'), f)
f = D + 'figS5_optimizer_convergence.xlsx'
for sh in ('0–20deg', '20–40deg', '40–60deg', '60–80deg', 'total EQE'):
    add('FigS6_' + sh.replace('–', '-').replace(' ', '_'), 'Fig. S6', sh, f'Best-so-far objective vs evaluation count, {sh}.', rows_of(f, sh), f)
f = 'fig3_efficiency_composition.xlsx'
add('FigS8a_points', 'Fig. S8', 'a', 'Total EQE vs 40-60 deg EQE: random and optimizer-visited designs.', rows_of(f, 'a_collapse_points'), f)
add('FigS8a_fit', 'Fig. S8', 'a', 'Linear fit (R^2 = 0.968).', rows_of(f, 'a_fit'), f)
add('FigS8b_dedicated', 'Fig. S8', 'b', 'Selectivity of band-dedicated optima vs natural composition.', rows_of(f, 'b_dedicated_vs_natural'), f)
add('FigS8c_points', 'Fig. S8', 'c', 'Selectivity points.', rows_of(f, 'c_selectivity_points'), f)
add('FigS8c_fits', 'Fig. S8', 'c', 'Selectivity fits.', rows_of(f, 'c_fits'), f)

# ---------------- order sheets: main figures, then supplementary figures, each by number
import re as _re
def _key(ws):
    if ws.title == 'README': return (-1, 0, '')
    mm = _re.match(r'Fig(S?)(\d+)', ws.title); return (1 if mm.group(1) else 0, int(mm.group(2)), ws.title)
wb._sheets.sort(key=_key)
index.sort(key=lambda r: (1 if r[1].startswith('Fig. S') else 0, int(_re.findall(r'\d+', r[1])[0]), r[0]))
# ---------------- README
readme['A1'] = 'Source data for all figures (main text Figs. 1-5, Supplement Figs. S1-S8)'; readme['A1'].font = Font(bold=True, size=12)
readme['A2'] = 'One sheet per panel. Angles are substrate (glass) polar angles in degrees; efficiencies are fractions of the power launched into the substrate unless noted; pp = percentage points.'
readme['A3'] = 'Structure classes: planar, scatter (g = asymmetry, t = reduced optical depth tau*), MLA hemi/free (benchmark optima per band), A_shape/B_aspect/C_jitter/C2_jitter (randomized lenses), W_slope (wavy surfaces).'
for j, t in enumerate(['sheet', 'figure', 'panel', 'content', 'generated from']):
    c = readme.cell(5, j + 1, t); c.font = HDR; c.fill = FILL
for i, row in enumerate(index):
    for j, v in enumerate(row):
        readme.cell(6 + i, j + 1, v).alignment = Alignment(wrap_text=(j == 3), vertical='top')
for col, w in zip('ABCDE', (26, 9, 7, 90, 42)):
    readme.column_dimensions[col].width = w
wb.save(OUT); print(OUT, len(index), 'sheets')
