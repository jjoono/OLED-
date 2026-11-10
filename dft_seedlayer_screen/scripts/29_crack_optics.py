"""Optical penalty of HATCN cracks: compare intact region (glass/HATCN30/Ag)
vs crack region (glass/Ag, HATCN missing). Area-weighted estimate of T loss.
"""
import numpy as np, openpyxl
from scipy.io import loadmat

wb = openpyxl.load_workbook(
    r"C:\Users\Junho\.claude\uploads\298d98bb-4f41-4a12-af8c-f1c865056809\80fe6f6b-Ag15nm_MgAg_8nm_260703.xlsx",
    data_only=True)
ws = wb["260703"]
lam=[];n2=[];k2=[]
for r in ws.iter_rows(min_row=5, values_only=True):
    if r[2] is None: continue
    lam.append(float(r[2])); n2.append(float(r[4])); k2.append(float(r[10]))
lam=np.array(lam); nAg=np.array(n2)+1j*np.array(k2)

m = loadmat(r"C:\Users\Junho\Dropbox\Linkstation\Simulation\LosslessEML_single_distribution\nk_JH_total.mat")
h = np.asarray(m["material"][0,0]["HATCN"]).ravel()
wl_h = np.arange(400, 400+len(h), 1.0)

def tmm(nl, dl, l):
    M=np.eye(2,dtype=complex)
    for j in range(len(nl)-1):
        a,b=nl[j],nl[j+1]; r=(a-b)/(a+b); t=2*a/(a+b)
        I=np.array([[1,r],[r,1]],dtype=complex)/t
        if j+1 < len(nl)-1:
            d=2*np.pi*nl[j+1]*dl[j+1]/l
            M=M@I@np.array([[np.exp(-1j*d),0],[0,np.exp(1j*d)]])
        else:
            M=M@I
    T=(np.real(nl[-1])/np.real(nl[0]))*abs(1/M[0,0])**2
    R=abs(M[1,0]/M[0,0])**2
    return T,R

print("HATCN 크랙의 광학 영향 — Ag 13.5 nm, 유리 기판, 수직입사")
print(f"{'nm':>4} | {'T_정상':>7} {'A_정상':>7} | {'T_크랙':>7} {'A_크랙':>7} | {'ΔT':>6}")
for p in (400,450,500,550,600,650,700):
    nA = np.interp(p,lam,nAg.real)+1j*np.interp(p,lam,nAg.imag)
    nH = (np.interp(p,wl_h,h.real,left=h[0].real,right=h[-1].real)
          +1j*np.interp(p,wl_h,h.imag,left=h[0].imag,right=h[-1].imag))
    Tn,Rn = tmm([1,nA,nH,1.52],[0,13.5,30,0],p)   # intact: HATCN present
    Tc,Rc = tmm([1,nA,1.52],[0,13.5,0],p)          # crack: HATCN absent
    print(f"{p:>4} | {100*Tn:7.1f} {100*(1-Tn-Rn):7.1f} | {100*Tc:7.1f} {100*(1-Tc-Rc):7.1f} | {100*(Tc-Tn):+6.1f}")

# area-weighted
print("\n크랙 면적률별 전체 T 변화 @550nm")
p=550
nA=np.interp(p,lam,nAg.real)+1j*np.interp(p,lam,nAg.imag)
nH=(np.interp(p,wl_h,h.real,left=h[0].real,right=h[-1].real)
    +1j*np.interp(p,wl_h,h.imag,left=h[0].imag,right=h[-1].imag))
Tn,_=tmm([1,nA,nH,1.52],[0,13.5,30,0],p)
Tc,_=tmm([1,nA,1.52],[0,13.5,0],p)
for f in (0.02,0.05,0.10,0.20):
    Teff=(1-f)*Tn+f*Tc
    print(f"  크랙 {100*f:4.0f}% : T = {100*Teff:.2f}%  (정상만이면 {100*Tn:.2f}%, Δ={100*(Teff-Tn):+.2f}%p)")
