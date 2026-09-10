"""Angle-resolved single-pass loss at the top electrode.
organic n=1.8 / metal d / capping n=2.1 / air, light incident from the organic
side.  In an OLED the emitter radiates into all internal angles, so what matters
is A(theta) for s and p separately - p is where the surface-plasmon channel is."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import numpy as np, json, csv, ellipsometry_fit as ef, ag_final as AG, tr_check as TC
TR=_os.path.join(ELLIPS_OUT, r'TR260820')
V3=json.load(open('ag_v3_result.json')); wl=np.array([550.0])
N_ORG=1.80; N_CAP=2.10; D_CAP=65.0

def RTA(N,d,th_org,pol):
    """theta measured inside the organic; returns (R,T,A) for one polarisation"""
    org=np.array([N_ORG+0j]); cap=np.array([N_CAP+0j]); air=np.array([1.0+0j])
    lay=[org,N,cap,air]; ds=[d,D_CAP]
    rp,rs=ef._tmm(wl,lay,ds,np.rad2deg(th_org))
    tp,ts=TC._t_amp(wl,lay,ds,np.rad2deg(th_org))
    s0=N_ORG*np.sin(th_org)
    qo=np.sqrt(N_ORG**2-s0**2+0j); qa=np.sqrt(1.0-s0**2+0j)
    if np.real(qa)<=0 or np.abs(np.imag(qa))>1e-9:      # beyond the escape cone
        Tt=0.0
    else:
        amp=(abs(tp[0])**2 if pol=='p' else abs(ts[0])**2)
        Tt=amp*np.real(qa)/np.real(qo) * (1.0 if pol=='s' else 1.0)
    R=(abs(rp[0])**2 if pol=='p' else abs(rs[0])**2)
    return R,Tt,max(0.0,1-R-Tt)

def load_nk(f,wt=550.):
    r=[x for x in csv.reader(open(f)) if x and not x[0].lstrip('\ufeff').startswith('#')]
    d=np.array([[float(v) for v in x[:3]] for x in r[1:] if x[0].replace('.','').replace('-','').isdigit()])
    return complex(np.interp(wt,d[:,0],d[:,1]),np.interp(wt,d[:,0],d[:,2]))

mc=np.genfromtxt(TR+r'\nk\Ag_McPeak.csv',delimiter=',',names=True); mw=mc[mc.dtype.names[0]]
if mw.max()<10: mw=mw*1000
N_bulk=np.array([complex(np.interp(550,mw,mc[mc.dtype.names[1]]),np.interp(550,mw,mc[mc.dtype.names[2]]))])
a=np.genfromtxt(TR+r'\nk\Ag8nm_on_HATCN5_measured.csv',delimiter=',',names=True)
N_tr=np.array([complex(np.interp(550,a['wavelength_nm'],a['n']),np.interp(550,a['wavelength_nm'],a['k']))])
N_se=AG.agN(np.array(V3['2-6']['p']),wl)
N_mg=np.array([load_nk('MgAg_nk.csv')])

crit=np.rad2deg(np.arcsin(1.0/N_ORG))
print('escape cone in the organic (n=1.8): %.1f deg'%crit)
print('\nabsorption in the metal, 8 nm, 550 nm, capping 65 nm')
print('%-20s | %s'%('n,k used','  '.join('%5.0f'%t for t in (0,20,40,50,60,70,80))))
for lab,N in [('Ag bulk (McPeak)',N_bulk),('Ag 8nm my SE',N_se),('Ag 8nm their T/R',N_tr),('Mg:Ag',N_mg)]:
    for pol in ('p','s'):
        row=[]
        for t in (0,20,40,50,60,70,80):
            R,T,A=RTA(N,8.0,np.deg2rad(t),pol); row.append('%5.1f'%(100*A))
        print('%-20s |  %s   (%s-pol)'%(lab if pol=='p' else '',' '.join(row),pol))
print()
print('same, but unpolarised average and for several thicknesses (their T/R n,k)')
print('%-8s | %s'%('d (nm)','  '.join('%5.0f'%t for t in (0,20,40,50,60,70,80))))
for d in (5,8,12,20):
    row=[]
    for t in (0,20,40,50,60,70,80):
        Rp,Tp,Ap=RTA(N_tr,float(d),np.deg2rad(t),'p')
        Rs_,Ts,As=RTA(N_tr,float(d),np.deg2rad(t),'s')
        row.append('%5.1f'%(100*(Ap+As)/2))
    print('%-8d |  %s'%(d,' '.join(row)))
