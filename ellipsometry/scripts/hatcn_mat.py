import re, numpy as np
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

D=str(DATA_DIR)
t=open(D+r'\HATCN_30.mat',encoding='utf-8-sig').read()

# split into the 3 B-Spline Fit Parms sections (uniaxial: x,y,z)
secs=re.split(r'start_B-Spline Fit Parms',t)[1:]
print('sections:',len(secs))

def parse_sec(s):
    einf=float(re.search(r"([\d.eE+-]+)\t[TF]\t-10\.0\t100\.0\tF\t'E Inf'",s).group(1))
    ira=re.search(r"([\d.eE+-]+)\t[TF]\t0\.0\t1000\.0\tF\t'IR Amp'",s)
    ir=float(ira.group(1)) if ira else 0.0
    nodes=re.findall(r"([\d.eE+-]+)\t[TF]\t[-\d.eE]+\t[\d.eE+-]+\t[TF]\t'spline_e2\(([\d.]+)\)'",s)
    E=np.array([float(b) for a,b in nodes]); V=np.array([float(a) for a,b in nodes])
    o=np.argsort(E)
    return E[o],V[o],einf,ir

# KK
Eg=np.linspace(0.4,9.0,1600); dE=Eg[1]-Eg[0]
Ej=Eg[None,:]; Ei=Eg[:,None]
with np.errstate(divide='ignore',invalid='ignore'):
    M=Ej/(Ej**2-Ei**2); S=1.0/(Ej**2-Ei**2)
np.fill_diagonal(M,0); np.fill_diagonal(S,0)
Ssum=S.sum(1); a_,b_=Eg[0],Eg[-1]
with np.errstate(divide='ignore'):
    Iana=(1/(2*Eg))*(np.log(np.abs((b_-Eg)/(b_+Eg)))-np.log(np.abs((a_-Eg)/(a_+Eg))))
Iana[~np.isfinite(Iana)]=0
def kk(e2):
    L=(e2+Eg*np.gradient(e2,Eg))/(2*Eg)
    return (2/np.pi)*(dE*(M@e2)-dE*Eg*e2*Ssum+dE*L+Eg*e2*Iana)
from scipy.interpolate import PchipInterpolator

wl=np.linspace(250,1000,400); E_wl=1239.84193/wl
out={}
for i,s in enumerate(secs):
    E,V,einf,ir=parse_sec(s)
    e2g=np.maximum(PchipInterpolator(E,V,extrapolate=False)(Eg),0.0); e2g[~np.isfinite(e2g)]=0
    e1g=einf+kk(e2g)
    e1=np.interp(E_wl,Eg,e1g); e2=np.interp(E_wl,Eg,e2g)
    N=np.sqrt((e1+1j*e2).astype(complex))
    out[i]=(N.real,N.imag)
    print('axis%d: einf=%.3f nodes=%d Erange %.3f-%.3f'%(i,einf,len(E),E[0],E[-1]))
    for tw in [400,450,550,633,800,950]:
        j=int(np.argmin(abs(wl-tw)))
        print('   %4dnm n=%.3f k=%.4f'%(wl[j],N[j].real,N[j].imag))

import csv
with open(str(OUT_DIR / 'HATCN_CE_nk.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['wl_nm','n_x','k_x','n_y','k_y','n_z','k_z'])
    for j in range(len(wl)):
        row=['%.2f'%wl[j]]
        for i in range(len(secs)):
            row+=['%.4f'%out[i][0][j],'%.4f'%out[i][1][j]]
        w.writerow(row)
print('saved HATCN_CE_nk.csv')
np.savez(str(OUT_DIR / 'HATCN_CE_nk.npz'),wl=wl,
         **{f'n{i}':out[i][0] for i in out},**{f'k{i}':out[i][1] for i in out})
