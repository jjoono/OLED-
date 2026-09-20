# -*- coding: utf-8 -*-
"""Raw data behind Fig. 3(b): the mode-resolved power budget against organic
thickness, for three substrate indices."""
import numpy as np, openpyxl, os
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

HERE = os.path.dirname(os.path.abspath(__file__))
BOLD = Font(bold=True)
F6, F1 = '0.000000', '0.0'

A = np.genfromtxt(os.path.join(HERE, 'fig3b_modes.csv'), delimiter=',')[1:]
G = np.genfromtxt(os.path.join(HERE, 'fig3b_ugrid_demo.csv'), delimiter=',')[1:]

wb = openpyxl.Workbook()
ws = wb.active
ws.title = 'README'
ws.column_dimensions['A'].width = 26
ws.column_dimensions['B'].width = 112
rows = [
 ('Fig. 3(b) - waveguide and SPP against organic thickness and substrate index', None),
 (None, None),
 ('Stack (bottom emitting, light leaves through the substrate)', None),
 ('  top / rear', 'air | Ag 100 nm   (McPeak, n = 0.044 + 3.819i at 550 nm)'),
 ('  organic', 'single isotropic layer, n = 1.80, thickness swept 10-500 nm'),
 ('  emitter', 'at the centre of the organic layer, PLQY = 1, isotropic orientation (horizontal fraction 2/3)'),
 ('  bottom electrode', 'ITO 50 nm, n = 1.86 + 0.003i'),
 ('  substrate', 'semi-infinite and incoherent, n = 1.50 / 1.65 / 1.80'),
 ('  wavelength', '550 nm, single wavelength'),
 (None, None),
 ('Method', None),
 ('  model', 'dipole transfer matrix (CPS): the dissipated power is resolved in the in-plane wavevector'),
 (None, 'u = k_x / (k0 n_organic) and partitioned into five channels that add up to 1 exactly'),
 ('  air', 'u < 1/n_organic and escaping the substrate, with the substrate treated incoherently'),
 ('  sub_confined', 'u < n_sub/n_organic, delivered into the substrate but trapped by total internal reflection at its front face'),
 ('  wg', 'n_sub/n_organic < u < 1: guided in the organic/ITO slab'),
 ('  spp', 'u > 1: surface plasmon polariton on the Ag, plus the rest of the evanescent near field'),
 ('  abs', 'absorbed inside the substrate cone (Ag + ITO)'),
 (None, None),
 ('  quadrature', 'the u integral is taken with u = sin(th) below the light line and u = sqrt(1+v^2) above it,'),
 (None, 'then composite Simpson with 12000 panels per sub-interval. The substitutions remove the'),
 (None, '1/sqrt(1-u^2) branch point at u = 1 analytically, and the channel boundaries sit exactly on'),
 (None, 'u = n_sub/n_organic and u = 1/n_organic instead of being snapped to a grid index.'),
 (None, 'Doubling the number of panels moves every entry by less than 1e-4.'),
 (None, None),
 ('  why this matters', 'the guided modes of the organic slab are poles of the integrand whose width in u is 1e-5..1e-3,'),
 (None, 'set only by the residual loss (Ag tail, ITO k = 0.003). A uniform grid with 1000 points per unit u'),
 (None, 'hits or misses them at random as the thickness is swept; see the sheet "u grid artefact", where'),
 (None, 'the waveguided channel scatters by up to 12 percentage points about the converged curve.'),
 (None, None),
 ('Sheets', None),
 ('  n_sub = 1.50 / 1.65 / 1.80', '2 nm steps, the curves plotted in the figure'),
 ('  disp_matrix (10 nm)', 'the same quantity on the 10:10:500 grid and in the column order of disp_matrix'),
 ('  u grid artefact', 'converged result against the uniform-grid result at N = 1000 and N = 3000'),
]
for i, (a, b) in enumerate(rows, 1):
    ws.cell(i, 1, a).font = BOLD if (a and not a.startswith(' ')) else Font()
    ws.cell(i, 2, b)
    ws.cell(i, 2).alignment = Alignment(wrap_text=False)

HDR = ['organic thickness (nm)', 'n_substrate', 'air', 'sub_confined', 'wg',
       'spp', 'abs', 'eta_sub = air + sub_confined', 'Purcell factor',
       'sum of the five channels']
for ns in (1.50, 1.65, 1.80):
    s = wb.create_sheet('n_sub = %.2f' % ns)
    for j, h in enumerate(HDR, 1):
        s.cell(1, j, h).font = BOLD
        s.column_dimensions[get_column_letter(j)].width = max(11, len(h) + 2)
    m = A[np.isclose(A[:, 1], ns)]
    for i, r in enumerate(m, 2):
        s.cell(i, 1, float(r[0])).number_format = F1
        s.cell(i, 2, float(r[1]))
        for j in range(3, 10):
            s.cell(i, j, float(r[j - 1])).number_format = F6
        s.cell(i, 10, '=SUM(C%d:G%d)' % (i, i)).number_format = F6
    s.freeze_panes = 'A2'

s = wb.create_sheet('disp_matrix (10 nm)')
H2 = ['d1 (nm)', 'n_substrate', 'EQE_air', 'EQE_sub_confined', 'EQE_wg',
      'EQE_spp', 'EQE_abs']
for j, h in enumerate(H2, 1):
    s.cell(1, j, h).font = BOLD
    s.column_dimensions[get_column_letter(j)].width = max(11, len(h) + 2)
i = 2
for ns in (1.50, 1.65, 1.80):
    m = A[np.isclose(A[:, 1], ns)]
    for r in m[np.isclose(np.mod(m[:, 0], 10.0), 0.0)]:
        s.cell(i, 1, float(r[0])).number_format = F1
        s.cell(i, 2, float(r[1]))
        for j in range(3, 8):
            s.cell(i, j, float(r[j - 1])).number_format = F6
        i += 1
s.freeze_panes = 'A2'

s = wb.create_sheet('u grid artefact')
H3 = ['organic thickness (nm)', 'wg  N=1000', 'spp  N=1000', 'sub  N=1000',
      'wg  N=3000', 'spp  N=3000', 'sub  N=3000',
      'wg  converged', 'spp  converged', 'sub  converged']
for j, h in enumerate(H3, 1):
    s.cell(1, j, h).font = BOLD
    s.column_dimensions[get_column_letter(j)].width = max(11, len(h) + 2)
for i, r in enumerate(G, 2):
    s.cell(i, 1, float(r[0])).number_format = F1
    for j in range(2, 11):
        s.cell(i, j, float(r[j - 1])).number_format = F6
s.freeze_panes = 'A2'

wb.save(os.path.join(HERE, 'fig3b_rawdata.xlsx'))
print('fig3b_rawdata.xlsx written: %d sweep rows' % len(A))
