"""J-V-L and L-EQE plots of representative devices.

Usage:  python3 plot_jvl_eqe.py bare=S2D2.csv freeform=S2D3.csv hemisphere=S2D4.csv [out=jvl_eqe.png]
Each input (.csv or .xlsx, first sheet) needs columns (header names, case-insensitive):
  V [V], J [mA/cm2], L [cd/m2], EQE [%]   (EQE = corrected EQE, i.e. 0-deg EQE x angular factor)
Optional column: EQE0 [%] (0-deg, Lambertian-assumed EQE) -> drawn dashed in panel (b).
Rows with L < 1 cd/m2 are kept for J-V but omitted from the EQE panel (below turn-on)."""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
C = {'bare': '#7d7c76', 'freeform': '#eb6834', 'hemisphere': '#2a78d6'}


def load(path):
    if path.lower().endswith(('.xlsx', '.xls')):
        from openpyxl import load_workbook
        rows = list(load_workbook(path, read_only=True, data_only=True).worksheets[0].iter_rows(values_only=True))
    else:
        import csv
        rows = list(csv.reader(open(path, encoding='utf-8-sig')))
    hdr = [str(h).strip().lower() for h in rows[0]]
    def col(*names):
        for n in names:
            for i, h in enumerate(hdr):
                if h.split(' ')[0].split('(')[0].split('[')[0] == n: return i
        return None
    idx = {k: col(*v) for k, v in dict(V=('v', 'voltage'), J=('j', 'current_density'), L=('l', 'luminance'), EQE=('eqe', 'corrected_eqe'), EQE0=('eqe0', 'front_eqe')).items()}
    data = {k: np.array([float(r[i]) for r in rows[1:] if r[i] not in (None, '')]) if i is not None else None for k, i in idx.items()}
    return data


args = dict(a.split('=', 1) for a in sys.argv[1:])
out = args.pop('out', 'jvl_eqe.png')
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw=dict(wspace=0.55))
a1b = a1.twinx()
for nm, path in args.items():
    d = load(path); c = C.get(nm, 'k')
    a1.semilogy(d['V'], np.abs(d['J']), '-', color=c, lw=1.2, label=nm)
    a1b.semilogy(d['V'], np.clip(d['L'], 1e-2, None), '--', color=c, lw=1.0)
    ok = d['L'] >= 1.0
    a2.semilogx(d['L'][ok], d['EQE'][ok], '-o', color=c, ms=2.5, lw=1.1, label=nm)
    if d.get('EQE0') is not None:
        a2.semilogx(d['L'][ok], d['EQE0'][ok], '--', color=c, lw=0.9)
a1.set_xlabel('voltage (V)'); a1.set_ylabel('current density (mA cm$^{-2}$), solid'); a1b.set_ylabel('luminance (cd m$^{-2}$), dashed')
a1.legend(fontsize=6, frameon=False, loc='upper left'); a1.set_title('(a) J–V–L', loc='left', fontsize=8)
a2.set_xlabel('luminance (cd m$^{-2}$)'); a2.set_ylabel('EQE (%)'); a2.set_ylim(bottom=0); a2.legend(fontsize=6, frameon=False, loc='lower left')
a2.set_title('(b) EQE–L (dashed: 0° Lambertian EQE)', loc='left', fontsize=8)
for a in (a1, a2): a.spines[['top']].set_visible(False)
a2.spines[['right']].set_visible(False)
fig.savefig(out, dpi=300, bbox_inches='tight'); print('wrote', out)
