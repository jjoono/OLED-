# -*- coding: utf-8 -*-
"""Raw data behind Fig. 2(j) alone: the spectrum-averaged round-trip loss against the angle
in the substrate for Al, Ag and the two dielectric stacks, with the DBR loss split into
absorbed and transmitted, the cos.sin weight, the layer-by-layer DBR thicknesses and the
flux- and spectrum-weighted numbers quoted in the text."""
import numpy as np, openpyxl, os, sys
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'sim', 'design_rule4'))
import fig2d as F

BOLD = Font(bold=True)
NSUB = 1.5
LAM = np.arange(430.0, 701.0, 1.0)
TH = np.radians(np.linspace(0.0, 89.9, 180))
SPEC = np.clip(np.interp(LAM, F.LAM, F.GREEN), 0, None); SPEC /= SPEC.sum()
CURVES = [('Al', 'Al'), ('Ag', 'Ag'), ('TCO + DBR', 'ITO 50 nm + chirped DBR'),
          ('TCO + DBR, plain', 'ITO 50 nm + plain quarter-wave DBR')]


def rt(name):
    S = F.stacks()[name]
    sel = np.isin(F.LAM, LAM)
    layers = [(m[sel], d) for m, d in S['layers']]
    return F.RT(layers, np.full(LAM.shape, NSUB, dtype=complex), S['exit'][sel], TH, LAM)


cols = {}
for name, lab in CURVES:
    R, T = rt(name)                      # (lam, th)
    A = np.clip(1.0 - R - T, 0.0, None)
    cols[lab] = (100 * (SPEC[:, None] * (1 - R)).sum(0), 100 * (SPEC[:, None] * A).sum(0),
                 100 * (SPEC[:, None] * T).sum(0))
w = np.cos(TH) * np.sin(TH); w /= w.max()

wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'README'
ws.column_dimensions['A'].width = 26; ws.column_dimensions['B'].width = 108
S_dbr = F.stacks()['TCO + DBR']; S_qw = F.stacks()['TCO + DBR, plain']
rows = [
 ('Fig. 2(j) — round-trip loss against the angle in the substrate', None),
 (None, 'Light arriving from the substrate on the OLED stack; 1 - R is what one round trip loses.'),
 (None, None),
 ('Stack', None),
 ('  common part', 'substrate (n = 1.50) / ITO 150 nm / 420 nm of non-absorbing organics (n = 1.8) / reflector.  The transparent electrode is fixed, only the reflector changes.'),
 ('  Al', 'measured n,k (project library), 100 nm'),
 ('  Ag', 'n,k from the project library (McPeak entry), 100 nm'),
 ('  chirped DBR', 'an ITO 50 nm cathode on the organics, then ten ZnS/LiF pairs whose thicknesses grow linearly through the stack by a factor 1.42 (ZnS 55 -> 77 nm, LiF 91 -> 127 nm); the three numbers minimise the flux- and spectrum-weighted loss over all substrate angles and the whole emission band.  Layer list on sheet "DBR layers".'),
 ('  plain DBR', 'the same ten pairs as a uniform quarter-wave stack at 550 nm (ZnS %.1f nm / LiF %.1f nm).  Dashed in (j); comparison only.' % (F.D_ZNS, F.D_LIF)),
 ('  ITO', 'n = 1.8636 + 0.0032285i at 550 nm — Koenig et al., ACS Nano 8, 6182 (2014), full dispersion (refractiveindex.info)'),
 (None, None),
 ('Quantities', None),
 ('  1 - R', 'the round-trip loss: absorbed plus, for the dielectric mirror, transmitted out of the back'),
 ('  absorbed / transmitted', 'the two parts of 1 - R.  T is zero for the metals; for the DBR it is non-zero only inside the escape cone to air, theta < asin(1/n_sub) = 41.8 deg, since beyond it the light is trapped by total internal reflection at the back surface.'),
 ('  spectrum averaging', 'every curve is averaged over the green emission spectrum (430-700 nm, 1 nm steps)'),
 ('  weight', 'cos(theta) sin(theta), the ergodic weight an angularly randomising outcoupling structure enforces, normalised to its maximum; drawn as the shaded area in (j)'),
 (None, None),
 ('Sheets', None),
 ('  j curves', 'the panel: angle in the substrate, the four 1 - R curves, the DBR split, the weight'),
 ('  DBR layers', 'layer-by-layer thicknesses of the chirped and the plain stacks'),
 ('  summary', 'the flux- and spectrum-weighted numbers quoted in the text'),
 (None, None),
 ('Provenance', 'figures/fig2_full/export_j.py on sim/design_rule4/fig2d.py (dispersive transfer matrix); repository jjoono/OLED-, branch claude/oled-efficiency-paper-ag-955zsk'),
]
for i, (a, b) in enumerate(rows, 1):
    ws.cell(i, 1, a).font = BOLD if (a and not a.startswith(' ')) else Font()
    ws.cell(i, 2, b)

s = wb.create_sheet('j curves')
hdr = ['angle in substrate (deg)', 'Al: 1-R (%)', 'Ag: 1-R (%)', 'ITO+chirped DBR: 1-R (%)',
       'ITO+chirped DBR: absorbed (%)', 'ITO+chirped DBR: transmitted (%)', 'ITO+plain DBR: 1-R (%)',
       'ITO+plain DBR: absorbed (%)', 'ITO+plain DBR: transmitted (%)', 'cos.sin weight (norm.)']
for j, h in enumerate(hdr, 1):
    s.cell(1, j, h).font = BOLD
    s.column_dimensions[get_column_letter(j)].width = 24 if j > 1 else 22
c_al, c_ag = cols['Al'][0], cols['Ag'][0]
c_d = cols['ITO 50 nm + chirped DBR']; c_q = cols['ITO 50 nm + plain quarter-wave DBR']
for i in range(len(TH)):
    vals = [np.degrees(TH[i]), c_al[i], c_ag[i], c_d[0][i], c_d[1][i], c_d[2][i], c_q[0][i], c_q[1][i], c_q[2][i], w[i]]
    for j, v in enumerate(vals, 1):
        s.cell(2 + i, j, float(v)).number_format = '0.000' if j > 1 else '0.00'
s.freeze_panes = 'A2'

s = wb.create_sheet('DBR layers')
for j, h in enumerate(['layer (from the substrate outward)', 'material', 'chirped stack (nm)', 'plain quarter-wave stack (nm)'], 1):
    s.cell(1, j, h).font = BOLD
d_c = [float(d) for _, d in S_dbr['layers']]
d_q = [float(d) for _, d in S_qw['layers']]
names = ['ITO (anode)', 'organics (n = 1.8, non-absorbing)', 'ITO (cathode)'] + \
        [('ZnS' if i % 2 == 0 else 'LiF') for i in range(len(d_c) - 3)]
for i, (nm, dc, dq) in enumerate(zip(names, d_c, d_q)):
    s.cell(2 + i, 1, i + 1); s.cell(2 + i, 2, nm)
    s.cell(2 + i, 3, dc).number_format = '0.0'; s.cell(2 + i, 4, dq).number_format = '0.0'
s.cell(3 + len(d_c), 1, 'beyond the last LiF: air.  Both stacks sit on the same substrate / ITO 150 nm / organics 420 nm / ITO 50 nm.')
for j in range(1, 5): s.column_dimensions[get_column_letter(j)].width = 34

s = wb.create_sheet('summary')
s['A1'] = 'Flux- and spectrum-weighted round-trip loss'; s['A1'].font = BOLD
s['A2'] = 'angles weighted cos.sin over 0-90 deg, wavelengths by the emission spectrum over 430-700 nm'
for j, h in enumerate(['reflector', 'absorbed (%)', 'transmitted (%)', 'total 1 − R (%)', 'substrate / emitter'], 1):
    s.cell(4, j, h).font = BOLD
r = 5
for ns, spec, tag in ((1.5, F.GREEN, 'n_sub 1.50, green emitter'), (1.8, F.ORANGE, 'n_sub 1.80, orange emitter')):
    for name, lab in CURVES:
        Wt = F.weighted(name, n_sub=ns, spectrum=spec, lam=np.arange(430.0, 701.0))
        s.cell(r, 1, lab)
        for j, v in enumerate((100 * Wt['absorbed'], 100 * Wt['transmitted'], 100 * Wt['loss']), 2):
            s.cell(r, j, float(v)).number_format = '0.00'
        s.cell(r, 5, tag); r += 1
    r += 1
s.column_dimensions['A'].width = 36
for j in range(2, 6): s.column_dimensions[get_column_letter(j)].width = 20
wb.save('fig2j_rawdata.xlsx')
print('saved fig2j_rawdata.xlsx; chirped layers', len(d_c), 'plain layers', len(d_q))
print('normal incidence: ' + ', '.join('%s %.2f' % (lab, cols[lab][0][0]) for _, lab in CURVES))
