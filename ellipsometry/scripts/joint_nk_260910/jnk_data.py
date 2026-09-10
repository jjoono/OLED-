"""Loaders for the 260819/260820 thin-Ag campaign: VASE on Si and absolute T/R on glass.

Two measurements of the same depositions on two substrates:

  SE   summary.xlsx, one sheet per Si piece, 5 angles 45-65 deg,
       675 points 245.8-1688.1 nm, blocks Psi/Delta | Re/Im(rho) | N/C/S |
       Intensity | % Depolarization.  Read from ELLIPS_DATA.
  T/R  the tracked 260820 campaign, dft_seedlayer_screen/data/TR_20260820/raw.
       Read from the raw exports and corrected with build_tra's own constant,
       so the numbers match the consolidated ALL_SAMPLES_TRA.csv wherever that
       file has them - and, unlike it, the two half-measured pieces keep the
       channel they do have (the delivered consolidation blanks MoOx5_bare and
       MoOx5_Ag12 entirely because each is missing its other channel).

Reflectance here is always the CORRECTED one.  The UMA collects only f = 0.837
of the substrate's back-surface beam at 6 deg, so the raw export is low by
(1 - f) Rg (T/100)^2 - 0.59 %p on bare glass, 0.2 %p at 12 nm of Ag.  Fitting
the raw R instead makes the substrate look absorbing and biases every
intensity-based thickness; scripts/tr_crosscheck/build_tra.py owns the
constant and this module imports it rather than repeating it.

The glass-piece and Si-piece numbers of one run are four apart (glass 1-2 <->
Si 1-6).  The glass side is fixed by the campaign handoff, which build_tra.py
carries; the Si side follows from the seed and Ag thicknesses in ag_load.SPEC,
and the T/R series agrees with it - absorption at 550 nm falls monotonically
14.9 -> 6.5 % across HATCN Ag 4 -> 12 nm and the Ag-free piece (1-9) reflects
like bare glass.

The workbook is not uniform: sheet 1-5 carries an extra leading column in the
depolarization block and 1-8 an extra header row, so the header row and every
block are located by name rather than by a fixed index.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), T/R csv
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import os, sys, numpy as np

_HERE = _os.path.dirname(_os.path.abspath(__file__))
sys.path.insert(0, _os.path.join(_HERE, '..', 'tr_crosscheck'))
import build_tra as BT                      # owns the raw layout and the R correction

XL    = _os.path.join(ELLIPS_DATA, 'summary.xlsx')
TRA   = _os.path.join(BT.DATA, 'ALL_SAMPLES_TRA.csv')
RAW   = BT.RAW
CACHE = _os.path.join(ELLIPS_OUT, 'jnk_data.npz')

ANG = np.array([45., 50., 55., 60., 65.])

# seed, nominal seed nm, nominal Ag nm, SE sheet (Si piece), T/R tag (glass piece),
# label in ALL_SAMPLES_TRA.csv
SAMPLES = [
    ('HATCN', 4,  0, '1-5',  None,   None),          # glass piece 1-1 was not measured
    ('HATCN', 4,  4, '1-6',  '1-2',  'HATCN5_Ag4'),
    ('HATCN', 4,  5, '1-7',  '1-3',  'HATCN5_Ag5'),
    ('HATCN', 4,  6, '1-8',  '1-4',  'HATCN5_Ag6'),
    ('HATCN', 4,  7, '2-5',  '2-1',  'HATCN5_Ag7'),
    ('HATCN', 4,  8, '2-6',  '2-2',  'HATCN5_Ag8'),
    ('HATCN', 4, 10, '2-7',  '2-3',  'HATCN5_Ag10'),
    ('HATCN', 4, 12, '2-8',  '2-4',  'HATCN5_Ag12'),
    ('MoOx',  5,  0, '1-13', '1-9',  'MoOx5_bare'),   # 1-9: R only
    ('MoOx',  5,  4, '1-14', '1-10', 'MoOx5_Ag4'),
    ('MoOx',  5,  5, '1-15', '1-11', 'MoOx5_Ag5'),
    ('MoOx',  5,  6, '1-16', '1-12', 'MoOx5_Ag6'),
    ('MoOx',  5,  7, '2-13', '2-9',  'MoOx5_Ag7'),
    ('MoOx',  5,  8, '2-14', '2-10', 'MoOx5_Ag8'),
    ('MoOx',  5, 10, '2-15', '2-11', 'MoOx5_Ag10'),
    ('MoOx',  5, 12, '2-16', '2-12', 'MoOx5_Ag12'),   # 2-12: T only
]
BY_SHEET = {s[3]: s for s in SAMPLES}
SHEETS   = [s[3] for s in SAMPLES]
SEEDS    = {'HATCN': '1-5', 'MoOx': '1-13'}          # the Ag-free reference of each seed


# ---------------------------------------------------------------- SE (xlsx)
def _find_blocks(rows):
    """Locate the header row and the first column of each named block."""
    for i, r in enumerate(rows[:10]):
        cols = {}
        for j, v in enumerate(r):
            if not isinstance(v, str):
                continue
            if v.startswith('Psi ('):              cols['psi'] = j
            elif v.startswith('Delta ('):          cols['del'] = j
            elif v.startswith('% Depolarization'): cols.setdefault('dep', j)
        if 'psi' in cols and 'del' in cols:
            return i, cols
    raise ValueError('no Psi/Delta header found')


def _read_xlsx():
    import openpyxl
    wb = openpyxl.load_workbook(XL, data_only=True, read_only=True)
    out = {}
    for sh in wb.sheetnames:
        rows = list(wb[sh].iter_rows(values_only=True))
        h, c = _find_blocks(rows)
        # Psi and Delta interleave from the Psi column: psi0 d0 psi1 d1 ...
        assert c['del'] == c['psi'] + 1, sh
        rec = []
        for r in rows[h + 1:]:
            try:
                wl = float(r[0])
            except (TypeError, ValueError, IndexError):
                continue
            psi = [float(r[c['psi'] + 2 * i]) for i in range(5)]
            dlt = [float(r[c['psi'] + 1 + 2 * i]) for i in range(5)]
            dep = []
            for i in range(5):
                try:
                    dep.append(float(r[c['dep'] + i]))
                except (TypeError, ValueError, IndexError, KeyError):
                    dep.append(np.nan)
            rec.append([wl] + psi + dlt + dep)
        out['SE_' + sh] = np.array(rec)
    return out


# ---------------------------------------------------------------- T/R (csv)
def _read_csv(path):
    wl, y = [], []
    for line in open(path, encoding='utf-8-sig'):
        p = line.strip().split(',')
        try:
            w, v = float(p[0]), float(p[1])
        except (ValueError, IndexError):
            continue
        wl.append(w); y.append(v)
    a = np.array(wl); b = np.array(y) / 100.0        # % -> fraction
    o = np.argsort(a)
    return np.column_stack([a[o], b[o]])


def _one_tr(tag):
    """Raw T and R of one piece, plus the corrected R, on the R grid.

    The bare-substrate transmittance export is named glass.csv, not glassT.csv.
    """
    pT = _os.path.join(RAW, 'glass.csv' if tag == 'glass' else '%sT.csv' % tag)
    pR = _os.path.join(RAW, '%sR.csv' % tag)
    T = _read_csv(pT) if os.path.exists(pT) else None
    Rm = _read_csv(pR) if os.path.exists(pR) else None
    out = {}
    if T is not None:
        out['T'] = T
    if Rm is not None:
        out['Rmeas'] = Rm
        if T is not None:
            Ti = np.interp(Rm[:, 0], T[:, 0], T[:, 1])
            out['R'] = np.column_stack([Rm[:, 0],
                                        Rm[:, 1] + BT.r_correction(Ti * 100.0) / 100.0])
    return out


def _read_tr():
    out = {}
    for tag in [s[4] for s in SAMPLES if s[4]] + ['glass']:
        for ch, arr in _one_tr(tag).items():
            out['TR_%s_%s' % (tag, ch)] = arr
    return out


def check_against_delivered(wl_lo=400.0, wl_hi=800.0):
    """Corrected R here vs the consolidated file, where that file has it."""
    A = np.genfromtxt(TRA, delimiter=',', names=True)
    wr = A['wavelength_nm']; o = np.argsort(wr)
    worst, worst_lab = 0.0, ''
    for smp in SAMPLES:
        tag, lab = smp[4], smp[5]
        key = '%s_Rcorr_pct' % lab if lab else None
        if not key or key not in A.dtype.names:
            continue
        ref = A[key][o]
        if not np.isfinite(ref).any():
            continue
        wl, mine = tr(tag, 'R', wl_lo, wl_hi)
        if wl is None:
            continue
        d = np.max(np.abs(100 * mine - np.interp(wl, wr[o], ref)))
        if d > worst:
            worst, worst_lab = d, lab
    return worst, worst_lab


# ---------------------------------------------------------------- public API
def _cache():
    if not os.path.exists(CACHE):
        d = _read_xlsx(); d.update(_read_tr())
        os.makedirs(_os.path.dirname(CACHE) or '.', exist_ok=True)
        np.savez_compressed(CACHE, **d)
    return np.load(CACHE)


def se(sheet, wl_min=260.0, wl_max=1080.0, step=1):
    """-> wl, Psi(nwl,5), Delta(nwl,5), depolarization(nwl,5), all in degrees/%

    step > 1 decimates the 675-point grid.  The spectra are smooth on the
    1.6-3 nm measurement spacing, so step=2 costs nothing but halves the fit
    time; it is what the oscillator fits use.
    """
    A = _cache()['SE_' + sheet]
    m = (A[:, 0] >= wl_min) & (A[:, 0] <= wl_max)
    A = A[m][::step]
    return A[:, 0], A[:, 1:6], A[:, 6:11], A[:, 11:16]


def tr(tag, ch, wl_min=400.0, wl_max=800.0):
    """-> wl, value (fraction).  Returns (None, None) if that channel is missing."""
    z = _cache(); key = 'TR_%s_%s' % (tag, ch)
    if key not in z.files:
        return None, None
    A = z[key]
    m = (A[:, 0] >= wl_min) & (A[:, 0] <= wl_max)
    return A[m, 0], A[m, 1]


def glass_tr(wl_min=400.0, wl_max=800.0, corrected=True):
    """Bare-substrate reference -> wl, T, R on the common grid."""
    z = _cache()
    T = z['TR_glass_T']; R = z['TR_glass_R' if corrected else 'TR_glass_Rmeas']
    wl = T[(T[:, 0] >= wl_min) & (T[:, 0] <= wl_max), 0]
    return wl, np.interp(wl, T[:, 0], T[:, 1]), np.interp(wl, R[:, 0], R[:, 1])


def ncs(psi_deg, del_deg):
    """N = cos2Psi, C = sin2Psi cosDelta, S = sin2Psi sinDelta."""
    p = np.deg2rad(psi_deg); d = np.deg2rad(del_deg)
    return np.cos(2 * p), np.sin(2 * p) * np.cos(d), np.sin(2 * p) * np.sin(d)


if __name__ == '__main__':
    wg, Tg, Rg = glass_tr()
    _, _, Rr = glass_tr(corrected=False)
    i = np.argmin(abs(wg - 550))
    print('bare glass @550 nm   T = %.2f %%   R = %.2f %% raw -> %.2f %% corrected'
          '   A = %+.2f %%p (raw %+.2f)'
          % (100 * Tg[i], 100 * Rr[i], 100 * Rg[i],
             100 * (1 - Tg[i] - Rg[i]), 100 * (1 - Tg[i] - Rr[i])))
    print('\n%-6s %-6s %4s %6s %6s | %5s %5s | %6s %6s %6s | %s'
          % ('sheet', 'seed', 'Ag', 'nSE', 'depol', 'tag', 'nTR',
             'T550', 'R550', 'A550', 'Psi@633 (45..65 deg)'))
    w, lab = check_against_delivered()
    print('corrected R here vs ALL_SAMPLES_TRA.csv: max %.4f %%p (%s)\n' % (w, lab))
    for seed, ds, da, sh, tag, lab in SAMPLES:
        wl, P, D, dep = se(sh)
        j = np.argmin(abs(wl - 633))
        wT, T = tr(tag, 'T') if tag else (None, None)
        wR, R = tr(tag, 'R') if tag else (None, None)
        g = lambda w, v: np.interp(550, w, v) if w is not None else np.nan
        t550, r550 = g(wT, T), g(wR, R)
        print('%-6s %-6s %4d %6d %6.2f | %5s %5s | %6s %6s %6s | %s'
              % (sh, seed, da, len(wl), np.nanmedian(dep),
                 tag or '-', len(wT) if wT is not None else '-',
                 '%6.2f' % (100 * t550) if wT is not None else '     -',
                 '%6.2f' % (100 * r550) if wR is not None else '     -',
                 '%6.2f' % (100 * (1 - t550 - r550)) if (wT is not None and wR is not None) else '     -',
                 ' '.join('%5.2f' % P[j, a] for a in range(5))))
