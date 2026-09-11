"""Build ALL_SAMPLES_TRA.csv from the raw Cary 6000i / UMA exports.

Input   data/TR_20260820/raw/{id}T.csv , {id}R.csv
        Each is a two-line header ("sampleNT" / "Wavelength (nm),%T,") followed by
        wavelength,value rows written from long to short wavelength.
Output  data/TR_20260820/ALL_SAMPLES_TRA.csv
        wavelength_nm, then per sample: T_pct, Rmeas_pct, Rcorr_pct, A_pct.

--------------------------------------------------------------------------------
The reflectance correction
--------------------------------------------------------------------------------
The UMA collects only part of the substrate's BACK-surface reflection: at 6 deg
that beam is displaced ~0.14 mm from the front-surface one and its edge misses
the aperture. Measured R is therefore low by the uncollected fraction of the
back-surface beam:

    R_meas = R_front + f * R_back        f = collected fraction
    R_corr = R_meas + (1 - f) * R_back

The back-surface beam has crossed the sample twice, so its weight scales with
the sample's own transmittance:

    R_back ~ Rg * (T/100)^2              Rg = glass/air reflectance at the back

which is why the correction is largest for bare glass (~0.65 %p) and shrinks as
the Ag gets thicker (0.44 %p at 5 nm, 0.20 %p at 12 nm) - exactly the trend the
handoff reports.

Calibration. The handoff measured f = 0.843 +- 0.024 on bare glass, validated by
the fact that applying it drives the bare-glass absorptance over 450-700 nm to
+0.01 +- 0.15 %p (i.e. to zero, correct for soda-lime in its transparent window)
while leaving the Fe3+/Fe2+ features intact.

CAVEAT ON PROVENANCE. The original consolidation script is not in the repository.
The constant below was recovered by least-squares against the delivered
ALL_SAMPLES_TRA.csv over 400-800 nm and reproduces its Rcorr column to
mean 0.004 %p / max 0.012 %p - far inside the sigma(A) = 0.15 %p measurement
noise. The implied collection factor is f = 0.837, inside the handoff's
0.843 +- 0.024. Run this file with --verify to re-check that agreement.

T needs no correction: at normal incidence the internally reflected components
stay collinear, so all of them are collected (bare glass measured T = 91.66 %
against a lossless-slab prediction of 91.47 %).
"""
import os as _os, sys, csv, numpy as np

HERE = _os.path.dirname(_os.path.abspath(__file__))
DATA = _os.environ.get('TR20260820_DIR',
                       _os.path.join(HERE, '..', '..', '..',
                                     'dft_seedlayer_screen', 'data', 'TR_20260820'))
RAW = _os.path.join(DATA, 'raw')
OUT = _os.path.join(DATA, 'ALL_SAMPLES_TRA.csv')

# sample id -> label, in the order the delivered CSV uses (handoff section 1)
SAMPLES = [('1-2', 'HATCN5_Ag4'), ('1-3', 'HATCN5_Ag5'), ('1-4', 'HATCN5_Ag6'),
           ('2-1', 'HATCN5_Ag7'), ('2-2', 'HATCN5_Ag8'), ('2-3', 'HATCN5_Ag10'),
           ('2-4', 'HATCN5_Ag12'),
           ('1-9', 'MoOx5_bare'), ('1-10', 'MoOx5_Ag4'), ('1-11', 'MoOx5_Ag5'),
           ('1-12', 'MoOx5_Ag6'), ('2-9', 'MoOx5_Ag7'), ('2-10', 'MoOx5_Ag8'),
           ('2-11', 'MoOx5_Ag10'), ('2-12', 'MoOx5_Ag12')]
# 1-9 has no T export and 2-12 has no R export; those columns stay blank.

N_GLASS = 1.5230                                   # soda-lime, visible
R_GLASS_BACK = ((N_GLASS - 1) / (N_GLASS + 1)) ** 2
F_COLLECTED = 0.8369                               # recovered; handoff quotes 0.843 +- 0.024
C_CORR = (1.0 - F_COLLECTED) * R_GLASS_BACK * 100.0  # %p per unit (T/100)^2


def read_raw(path):
    """two header lines, then wavelength,value; returns arrays sorted ascending"""
    wl, v = [], []
    with open(path, encoding='utf-8-sig') as fh:
        for i, line in enumerate(fh):
            if i < 2:
                continue
            p = line.split(',')
            try:
                wl.append(float(p[0])); v.append(float(p[1]))
            except (ValueError, IndexError):
                continue
    wl = np.array(wl); v = np.array(v); o = np.argsort(wl)
    return wl[o], v[o]


def r_correction(T_pct):
    """percentage points to ADD to the measured reflectance"""
    return C_CORR * (T_pct / 100.0) ** 2


def build(out_path=None):
    grid = None
    cols = {}
    for sid, lab in SAMPLES:
        for kind in ('T', 'R'):
            p = _os.path.join(RAW, '%s%s.csv' % (sid, kind))
            if not _os.path.exists(p):
                print('  missing export: %s%s.csv  (%s)' % (sid, kind, lab))
                continue
            wl, v = read_raw(p)
            if grid is None:
                grid = wl
            cols[(lab, kind)] = np.interp(grid, wl, v, left=np.nan, right=np.nan)

    rows = {'wavelength_nm': grid}
    header = ['wavelength_nm']
    for _, lab in SAMPLES:
        T = cols.get((lab, 'T')); R = cols.get((lab, 'R'))
        nan = np.full(len(grid), np.nan)
        T = nan if T is None else T
        R = nan if R is None else R
        Rc = R + r_correction(T)
        A = 100.0 - T - Rc
        for suf, arr in (('T_pct', T), ('Rmeas_pct', R), ('Rcorr_pct', Rc), ('A_pct', A)):
            header.append('%s_%s' % (lab, suf)); rows['%s_%s' % (lab, suf)] = arr

    order = np.argsort(-grid)                       # delivered file runs long -> short
    dest = out_path or OUT
    with open(dest, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh); w.writerow(header)
        for i in order:
            w.writerow([('' if not np.isfinite(rows[h][i]) else
                         ('%g' % rows[h][i] if h == 'wavelength_nm' else '%.4f' % rows[h][i]))
                        for h in header])
    print('wrote %s  (%d wavelengths, %d samples)' % (dest, len(grid), len(SAMPLES)))
    return rows, grid


def verify(rows, grid):
    """compare against the delivered file, if it is still on disk under another name"""
    ref = _os.path.join(DATA, 'ALL_SAMPLES_TRA.csv')
    D = np.genfromtxt(ref, delimiter=',', names=True)
    wr = D['wavelength_nm']; o = np.argsort(wr)
    worst = 0.0; worst_lab = ''
    for _, lab in SAMPLES:
        k = lab + '_Rcorr_pct'
        if k not in D.dtype.names:
            continue
        a = np.interp(grid, wr[o], D[k][o])
        b = rows[k]
        m = np.isfinite(a) & np.isfinite(b) & (grid >= 400) & (grid <= 800)
        if m.sum() == 0:
            continue
        e = np.max(np.abs(a[m] - b[m]))
        if e > worst:
            worst, worst_lab = e, lab
    print('max |Rcorr rebuilt - delivered| = %.4f %%p  (%s)' % (worst, worst_lab))


if __name__ == '__main__':
    # --verify rebuilds to a scratch file so the delivered CSV is never clobbered
    if '--verify' in sys.argv:
        rows, grid = build(OUT + '.rebuilt')
        verify(rows, grid)
        _os.remove(OUT + '.rebuilt')
    else:
        rows, grid = build()
