"""Is the measured Intensity channel reproduced by the model?
CompleteEASE's reflection Intensity ~ (Rp+Rs)/2 times an UNKNOWN instrument
throughput factor, so only the spectral SHAPE is meaningful unless a calibrated
reflectance reference was taken."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import numpy as np, openpyxl, json
import ce_fit as cf, ce_fit2 as f2, ellipsometry_fit as ef

XL = _os.path.join(ELLIPS_DATA, r'se추출.xlsx')

def load_int(sheet):
    wb = openpyxl.load_workbook(XL, data_only=True, read_only=True)
    ws = wb[sheet]; rows = []
    for r in ws.iter_rows(min_row=4, values_only=True):
        try: rows.append([float(r[0]), float(r[28]), float(r[29]), float(r[30]),
                          float(r[33]), float(r[34]), float(r[35])])
        except (TypeError, ValueError, IndexError): pass
    A = np.array(rows)
    return A[:,0], A[:,1:4], A[:,4:7]

def model_R(p, wl):
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    Nf = f2.film_N(p, wl); Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf)); amb = np.ones_like(Nf)
    out = []
    for an in cf.ANG0 + p[2]:
        rp, rs = ef._tmm(wl, [amb, Nr, Nf, ox, si], [p[1], p[0], cf.D_OX], an)
        out.append((np.abs(rp)**2 + np.abs(rs)**2) / 2.0)
    return np.array(out).T

R = json.load(open('ce_final_result.json'))
for sh in ['#1', '#3']:
    wl, I, dep = load_int(sh)
    p = np.array(R[sh]['p'])
    Rm = model_R(p, wl)
    m = (wl >= 340) & (wl <= 1080)
    print('=== %s ===' % sh)
    print('  wl      I_meas(65/70/75)          R_model(65/70/75)        ratio meas/model      depol%%')
    for t in [250, 350, 450, 550, 700, 900, 1080, 1400, 1650]:
        i = np.argmin(abs(wl-t))
        rat = I[i]/Rm[i]
        print('  %6.0f  %5.3f %5.3f %5.3f   %5.3f %5.3f %5.3f   %5.2f %5.2f %5.2f   %4.1f %4.1f %4.1f'
              % (wl[i], *I[i], *Rm[i], *rat, *dep[i]))
    rat = (I/Rm)[m]
    print('  ratio over 340-1080 nm: mean %.3f, spread %.1f%% (65deg %.3f+-%.3f)'
          % (rat.mean(), 100*rat.std()/rat.mean(), rat[:,0].mean(), rat[:,0].std()))
    # shape-only agreement: scale each angle by its own best constant
    for a in range(3):
        s = (I[m,a]*Rm[m,a]).sum()/(Rm[m,a]**2).sum()
        rms = np.sqrt(np.mean((I[m,a]-s*Rm[m,a])**2))/I[m,a].mean()
        print('    %2.0fdeg  best scale %.3f  ->  shape mismatch %.1f%% rms' % (cf.ANG0[a], s, 100*rms))
