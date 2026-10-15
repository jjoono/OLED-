"""Regenerate the delivered n,k CSV directly from the same parameters that the
.mod files carry, so CSV and .mod describe identical optical constants."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import json, csv, numpy as np
R=json.load(open(_os.path.join(ELLIPS_OUT, r'genosc_params.json')))
LAB={'1':'VendorA_ITO_1','2':'VendorA_ITO_2','3':'VendorB_IZO_1','4':'VendorB_IZO_2','5':'VendorA_ITO_2pctO2'}

EG=np.linspace(0.30,14.0,2400); dE=EG[1]-EG[0]
Ej=EG[None,:]; Ei=EG[:,None]
with np.errstate(divide='ignore',invalid='ignore'):
    MK=Ej/(Ej**2-Ei**2); SK=1.0/(Ej**2-Ei**2)
np.fill_diagonal(MK,0.0); np.fill_diagonal(SK,0.0)
Ssum=SK.sum(1); a_,b_=EG[0],EG[-1]
with np.errstate(divide='ignore'):
    Iana=(1/(2*EG))*(np.log(np.abs((b_-EG)/(b_+EG)))-np.log(np.abs((a_-EG)/(a_+EG))))
Iana[~np.isfinite(Iana)]=0.0
def kk(e2):
    L=(e2+EG*np.gradient(e2,EG))/(2*EG)
    return (2/np.pi)*(dE*(MK@e2)-dE*EG*e2*Ssum+dE*L+EG*e2*Iana)
def tl_e2(A,E0,C,Eg):
    den=(EG**2-E0**2)**2+C**2*EG**2
    return np.where(EG>Eg, A*E0*C*(EG-Eg)**2/(den*np.maximum(EG,1e-9)), 0.0)
def gau_e2(A,Ec,sig):
    return A*(np.exp(-((EG-Ec)/sig)**2)-np.exp(-((EG+Ec)/sig)**2))
def nk(p,E):
    e2=np.maximum(tl_e2(p[2],p[3],p[4],p[5])+gau_e2(p[6],p[7],p[8]),0.0)
    e1=p[0]+kk(e2)-p[1]/(EG**2+0.10**2)
    e2=e2+p[1]*0.10/(EG*(EG**2+0.10**2))
    N=np.sqrt((np.interp(E,EG,e1)+1j*np.interp(E,EG,np.maximum(e2,0))).astype(complex))
    N=np.where(N.imag<0,np.conj(N),N)
    return N.real,np.maximum(N.imag,0)

wl=np.arange(300.0,1101.0,2.0); E=1239.84193/wl
cols={}
for s in '12345':
    n,k=nk(R[s]['p'],E); cols[s]=(n,k)

notes=[
 ['# Optical constants (n,k) - sputtered ITO / IZO on Si + native oxide, batch 260813  [v3]'],
 ['# Generated from the SAME Gen-Osc parameters as the accompanying .mod files.'],
 ['# Model: Einf + Drude + Tauc-Lorentz + Gaussian, Kramers-Kronig consistent, eps2 >= 0.'],
 ['#'],
 ['# Geometry (thickness set by the independent glass-baseline transmittance, not SE alone;'],
 ['#  SE self-consistency is flat over d = 48-53 nm because of the n*d alias):'],
 ['#   VendorA_ITO_1/2 = 51 nm (rough 3 nm)   VendorB_IZO_1/2 = 42 nm (rough 2 nm)'],
 ['#   VendorA_ITO_2pctO2 = 42 nm (rough 2 nm)     native oxide 3 nm, Si substrate'],
 ['#'],
 ['# Uncertainty on n: +-0.05  (thickness alias +-0.03, substrate model +-0.03)'],
 ['# Reliable range 360-1080 nm. Outside it the values are model extrapolation:'],
 ['#   the UV band edge sits at the edge of the measured range, so E0_TL / C_TL are weakly determined.'],
 ['#'],
 ['# RELIABILITY'],
 ['#   VendorA_ITO_1/2    : n OK, k OK   (the pair agrees once the thickness is set correctly)'],
 ['#   VendorB_IZO_1/2   : n OK, k OK   (k < 0.01, at the detection limit - not literally zero)'],
 ['#   VendorA_ITO_2pctO2 : DO NOT USE   (alignment fault, needs 0.8 deg angle offset; T mismatch 8.7 %p)'],
 ['#'],
]
out=_os.path.join(ELLIPS_OUT, r'ITO_IZO_final_nk.csv')
with open(out,'w',newline='') as f:
    w=csv.writer(f)
    for nline in notes: w.writerow(nline)
    hdr=['wl_nm']
    for s in '12345': hdr+=['n_'+LAB[s],'k_'+LAB[s]]
    w.writerow(hdr)
    for i in range(len(wl)):
        row=['%.1f'%wl[i]]
        for s in '12345':
            row+=['%.4f'%cols[s][0][i],'%.4f'%cols[s][1][i]]
        w.writerow(row)
print('rewrote %s : %d rows, %.0f-%.0f nm'%(out,len(wl),wl[0],wl[-1]))
for tw in [450,550,633,800]:
    i=int(np.argmin(np.abs(wl-tw)))
    print(' %4dnm  '%tw+'  '.join('%s n=%.3f k=%.4f'%(LAB[s][:11],cols[s][0][i],cols[s][1][i]) for s in '12345'))
