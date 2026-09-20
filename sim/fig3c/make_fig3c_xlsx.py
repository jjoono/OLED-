# -*- coding: utf-8 -*-
"""Raw data behind Fig. 3(c): the ETL's out-of-plane index against the plasmon.

Sheets ci and cii are exactly what the two panels plot; the rest is supporting.
"""
import numpy as np, openpyxl, os
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter
import nspp, materials as M

HERE = os.path.dirname(os.path.abspath(__file__))
BOLD = Font(bold=True)
NE = [1.80, 1.70, 1.60, 1.50]
TARGETS = [0.40, 0.30, 0.20, 0.15, 0.10]
D_FIX = [40, 60, 80, 100, 150]

SP = np.genfromtxt(os.path.join(HERE, 'fig3c_spectrum.csv'), delimiter=',',
                   names=True, dtype=None, encoding='utf-8')
SW = np.genfromtxt(os.path.join(HERE, 'fig3c_sweep.csv'), delimiter=',',
                   names=True, dtype=None, encoding='utf-8')
FAM = SW[SW['series'] == 'no1.80']

wb = openpyxl.Workbook()
ws = wb.active
ws.title = 'README'
ws.column_dimensions['A'].width = 26
ws.column_dimensions['B'].width = 110
rows = [
 ('Fig. 3(c) - a low out-of-plane ETL index buys a thinner layer', None),
 (None, None),
 ('Stack', None),
 ('  top / rear', 'air | Ag 100 nm  (McPeak, 0.04382 + 3.81898i at 550 nm)'),
 ('  ETL', 'uniaxial, n_o = 1.80 fixed, n_e = 1.80 / 1.70 / 1.60 / 1.50, thickness swept'),
 ('  EML', '1.80 isotropic, 20 nm, dipole at the centre, PLQY = 1, isotropic orientation'),
 ('  HTL', '1.80 isotropic, 50 nm'),
 ('  bottom electrode', 'ITO 50 nm, 1.86362 + 0.00323i  (Koenig 2014)'),
 ('  substrate', 'semi-infinite and incoherent, n = 1.80'),
 ('  wavelength', '550 nm'),
 (None, None),
 ('Method', None),
 ('  model', 'dipole transfer matrix (CPS) with uniaxial layers; five channels summing to 1:'),
 (None, 'air + substrate-confined + waveguided + SPP/evanescent + absorbed'),
 ('  quadrature', 'u = sin(th) below the light line, u = sqrt(1+v^2) above it, composite Simpson'),
 (None, 'with 12000 panels per sub-interval. Reproduces the Fig. 3(b) solver to 3e-16 in the'),
 (None, "isotropic limit and the author's TMF_birefringence_whole.m to 6e-9 for a uniaxial stack."),
 ('  n_sub = n_EML = 1.80', 'the substrate and organic light lines coincide, so there is no waveguided channel:'),
 (None, 'power below k_x/k0 = 1.80 reaches the substrate, power above it does not.'),
 (None, None),
 ('What the panels show', None),
 ('  ci', 'TM dissipation density per unit k_x/k0 at d_ETL = 60 nm. The abscissa is the in-plane'),
 (None, 'wavevector in units of k0, so a peak sits at its mode effective index. The plasmon is at'),
 (None, '2.05 / 1.96 / 1.88 / 1.81 for n_e = 1.80 / 1.70 / 1.60 / 1.50: a lower out-of-plane'),
 (None, 'index moves it towards the substrate light line (1.80) and weakens it.'),
 ('  cii', 'the power left in the SPP against ETL thickness. 20 % is reached at'),
 (None, '110 / 104 / 90 / 59 nm, i.e. the same suppression with a thinner transport layer.'),
 (None, None),
 ('Caveat', 'the thicknesses depend on the EML and HTL indices and on wavelength; the trend'),
 (None, 'does not. The analytic plasmon index at a metal / uniaxial interface is'),
 (None, '(k_SPP/k0)^2 = eps_m (eps_o - eps_m)/(eps_o - eps_m^2/eps_e), giving 2.041 / 1.922 /'),
 (None, '1.804 / 1.687 in the thick-ETL limit for the four n_e above.'),
 (None, None),
 ('Sheets', None),
 ('  ci spectrum', 'TM and TE density vs k_x/k0, at d_ETL = 60, 100 and 150 nm'),
 ('  cii SPP vs thickness', 'the five channels vs d_ETL for the four n_e'),
 ('  summary', 'thickness at a given SPP target, and the budget at fixed thicknesses'),
 ('  full n_e sweep', 'the same, for n_e = 1.40 to 1.80 in steps of 0.02'),
 ('  measured ETLs', 'B3PyMPM / B4PyMPM / TPBi / TCTA with their own n_o (not plotted)'),
]
for i, (a, b) in enumerate(rows, 1):
    ws.cell(i, 1, a).font = BOLD if (a and not a.startswith(' ')) else Font()
    ws.cell(i, 2, b)


def write(name, arr, fmt=None, cols=None):
    s = wb.create_sheet(name)
    names = cols or arr.dtype.names
    for j, h in enumerate(names, 1):
        s.cell(1, j, h).font = BOLD
        s.column_dimensions[get_column_letter(j)].width = max(11, len(h) + 2)
    for i, r in enumerate(arr, 2):
        for j, h in enumerate(names, 1):
            v = r[h]
            c = s.cell(i, j, v.item() if hasattr(v, 'item') else v)
            if not isinstance(c.value, str):
                c.number_format = (fmt or {}).get(h, '0.000000')
    s.freeze_panes = 'A2'
    return len(arr)


keys = ['n_e=%.2f' % x for x in NE]
sel = np.isin(SP['series'], keys)
n1 = write('ci spectrum', SP[sel],
           {'d_ETL_nm': '0', 'n_eff': '0.00000', 'TM': '0.000E+00',
            'TE': '0.000E+00'})
sel = np.isin(np.round(FAM['n_e'], 2), NE)
n2 = write('cii SPP vs thickness', FAM[sel], {'d_ETL_nm': '0.0'})

s = wb.create_sheet('summary')
s.cell(1, 1, 'ETL thickness (nm) at which the SPP channel falls to a given level').font = BOLD
s.cell(2, 1, 'n_e').font = BOLD
s.cell(2, 2, 'n_SPP (thick-ETL limit)').font = BOLD
for j, t in enumerate(TARGETS, 3):
    s.cell(2, j, 'SPP = %.0f %%' % (100 * t)).font = BOLD
for i, ne in enumerate(NE, 3):
    m = FAM[np.isclose(FAM['n_e'], ne)]
    o = np.argsort(m['d_ETL_nm'])
    d, sp = m['d_ETL_nm'][o], m['spp'][o]
    s.cell(i, 1, ne).number_format = '0.00'
    s.cell(i, 2, float(nspp.n_spp(M.AG, 1.8, ne))).number_format = '0.0000'
    for j, t in enumerate(TARGETS, 3):
        k = np.where(sp <= t)[0]
        v = (np.interp(t, [sp[k[0]], sp[k[0] - 1]], [d[k[0]], d[k[0] - 1]])
             if len(k) and k[0] else float('nan'))
        s.cell(i, j, float(v)).number_format = '0.0'
r0 = 3 + len(NE) + 2
s.cell(r0, 1, 'SPP fraction / eta_sub at fixed ETL thickness').font = BOLD
s.cell(r0 + 1, 1, 'n_e').font = BOLD
for j, d in enumerate(D_FIX, 2):
    s.cell(r0 + 1, j, 'SPP @ %d nm' % d).font = BOLD
    s.cell(r0 + 1, j + len(D_FIX), 'eta_sub @ %d nm' % d).font = BOLD
for i, ne in enumerate(NE, r0 + 2):
    m = FAM[np.isclose(FAM['n_e'], ne)]
    o = np.argsort(m['d_ETL_nm'])
    s.cell(i, 1, ne).number_format = '0.00'
    for j, d in enumerate(D_FIX, 2):
        s.cell(i, j, float(np.interp(d, m['d_ETL_nm'][o], m['spp'][o]))
               ).number_format = '0.0000'
        s.cell(i, j + len(D_FIX), float(np.interp(d, m['d_ETL_nm'][o],
               m['eta_sub'][o]))).number_format = '0.0000'
for j in range(1, 2 * len(D_FIX) + 2):
    s.column_dimensions[get_column_letter(j)].width = 15

n3 = write('full n_e sweep', FAM, {'d_ETL_nm': '0.0'})
sel = np.isin(SW['series'], list(M.ETL))
n4 = write('measured ETLs', SW[sel], {'d_ETL_nm': '0.0'})

wb.save(os.path.join(HERE, 'fig3c_rawdata.xlsx'))
print('fig3c_rawdata.xlsx : ci %d, cii %d, sweep %d, measured %d rows'
      % (n1, n2, n3, n4))
