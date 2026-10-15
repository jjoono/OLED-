import numpy as np, ellipsometry_fit as ef
from scipy.optimize import minimize
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

base=str(DATA_DIR)
fp=base+r'\#2.xlsx'; depf=base+r'\#2_depol.xlsx'
angles=[65,70,75]
Ng=ef.load_grahil_nk(str(OUT_DIR / 'GraHIL_genosc_nk.csv'))
wl,psi_m,del_m=ef.load_data_xlsx(fp,angles,380,1688)
Wd=ef.depol_weights(wl,angles,depf,depol_cut=3.0,depol_soft=1.5)

# Graded surface: M sublayers, void fraction linear from f_bot(bottom) to f_top(top)
def calc_graded(wlq, Np, d_bulk, d_grade, f_top, f_bot=0.0, M=8):
    n_air=np.ones(len(wlq),dtype=complex)
    layers=[n_air]; d_list=[]
    # top (high void) first
    for m in range(M):
        frac_void = f_top + (f_bot-f_top)*(m+0.5)/M  # decreasing void with depth
        Nsub=ef.bruggeman_ema(Np, n_air, 1.0-frac_void)  # perovskite fraction = 1-void
        layers.append(Nsub); d_list.append(d_grade/M)
    layers.append(Np); d_list.append(d_bulk)
    layers.append(ef.n_SiO2(wlq)); d_list.append(2.0)
    layers.append(ef.n_Si(wlq))
    po,do=[],[]
    for ang in angles:
        rp,rs=ef._tmm(wlq,layers,d_list,ang); rho=rp/rs
        po.append(np.degrees(np.arctan(np.abs(rho)))); do.append(np.degrees(np.angle(rho)))
    return np.array(po).T, np.array(do).T

# per-wl free (n,k) residual at blue wavelengths for various graded params
blue_idx=[int(np.argmin(np.abs(wl-t))) for t in [430,460,500,540]]
def blue_res(d_bulk,d_grade,f_top):
    tot=0
    for i in blue_idx:
        w=wl[i:i+1]
        def f(x):
            n,k=x
            if n<1.3 or n>2.9 or k<0 or k>1.5: return 1e3
            Np=np.array([n+1j*k])
            pc,dc=calc_graded(w,Np,d_bulk,d_grade,f_top)
            dd=(dc[0]-del_m[i]+180)%360-180
            return np.sqrt(np.mean((pc[0]-psi_m[i])**2+dd**2))
        best=1e9
        for n0 in [1.8,2.2,2.6]:
            for k0 in [0.2,0.5,0.8]:
                r=minimize(f,[n0,k0],method='Nelder-Mead',options={'xatol':1e-4,'fatol':1e-6})
                if r.fun<best: best=r.fun
        tot+=best
    return tot/len(blue_idx)

print('Reference: single 80nm EMA blue per-wl res ~ 11 (from earlier)')
print('Graded surface scan (blue 430-540nm mean per-wl free-(n,k) res):')
for dg in [80,120,160]:
    for ft in [0.5,0.7,0.9]:
        # keep total ~270: bulk = 270 - dg*avg_void... approx bulk=270-dg*0.5*(ft/2)
        db=270-0.4*dg
        print('  d_grade=%3d f_top=%.1f bulk=%.0f: blue_res=%.2f'%(dg,ft,db,blue_res(db,dg,ft)))
