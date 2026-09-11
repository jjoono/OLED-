"""Loader for the 260819 ThinAg summary workbook.
16 samples, 5 angles (45-65 deg), 675 points 245.8-1688.1 nm, with % depolarization."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import numpy as np, openpyxl, os

XL = _os.path.join(ELLIPS_DATA, r'summary.xlsx')
CACHE = _os.path.join(ELLIPS_OUT, r'ag260819.npz')
ANG = np.array([45., 50., 55., 60., 65.])
SPEC = {  # sheet: (seed, seed nominal nm from the quartz monitor, Ag nominal nm)
          # Both seeds were deposited to a QCM reading of 5 nm; the fits require
          # 6.4-7.5 nm (HATCN) and 6.9-8.0 nm (MoOx). See README_260820_session.md.
 '1-5':('HATCN',5,0), '1-6':('HATCN',5,4), '1-7':('HATCN',5,5), '1-8':('HATCN',5,6),
 '2-5':('HATCN',5,7), '2-6':('HATCN',5,8), '2-7':('HATCN',5,10),'2-8':('HATCN',5,12),
 '1-13':('MoOx',5,0), '1-14':('MoOx',5,4), '1-15':('MoOx',5,5), '1-16':('MoOx',5,6),
 '2-13':('MoOx',5,7),'2-14':('MoOx',5,8),'2-15':('MoOx',5,10),'2-16':('MoOx',5,12)}
ORDER = ['1-5','1-6','1-7','1-8','2-5','2-6','2-7','2-8',
         '1-13','1-14','1-15','1-16','2-13','2-14','2-15','2-16']

def _read():
    wb = openpyxl.load_workbook(XL, data_only=True, read_only=True)
    D = {}
    for s in wb.sheetnames:
        rows = []
        for r in wb[s].iter_rows(min_row=4, values_only=True):
            try:
                wl = float(r[0])
            except (TypeError, ValueError):
                continue
            pd_ = [float(r[1+2*i]) for i in range(5)] + [float(r[2+2*i]) for i in range(5)]
            dep = []
            for j in range(50, 55):
                try: dep.append(float(r[j]))
                except (TypeError, ValueError, IndexError): dep.append(np.nan)
            rows.append([wl]+pd_+dep)
        A = np.array(rows)
        D[s] = A
    return D

def load(sheet):
    if not os.path.exists(CACHE):
        D = _read()
        np.savez_compressed(CACHE, **D)
    z = np.load(CACHE)
    A = z[sheet]
    return A[:,0], A[:,1:6], A[:,6:11], A[:,11:16]     # wl, Psi(5), Delta(5), depol(5)

if __name__ == '__main__':
    print('%-6s %-6s %4s %4s   %s' % ('sheet','seed','seed','Ag','depolarization %% (median over 300-1080 nm), per angle'))
    for s in ORDER:
        wl, P, Dl, dep = load(s)
        m = (wl>=300)&(wl<=1080)
        seed, ds, da = SPEC[s]
        print('%-6s %-6s %4d %4d   %s' % (s, seed, ds, da,
              '  '.join('%5.2f'%np.nanmedian(dep[m,i]) for i in range(5))))
