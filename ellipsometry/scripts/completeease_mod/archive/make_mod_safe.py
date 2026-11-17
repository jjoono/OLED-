"""Safe .mod generator.

The previous files added a 3rd oscillator (Gaussian); CompleteEASE evidently
refused them.  Here the template's structure is preserved EXACTLY - same
oscillator count (2: Drude + Tauc-Lorentz), same number of lines, same tabs -
and only numeric fields are substituted.

Because the Gaussian is dropped, the dispersion is re-fitted with Drude + TL
only, against the same per-wavelength (chain) dispersion.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, os, math, json, csv
import numpy as np
from scipy.optimize import least_squares

HBAR_EVS=6.582119569e-16; HBAR_EVFS=0.6582119569; EPS0=8.8541878128e-12
TEMPLATE=_os.path.join(ELLIPS_DATA, r'VendorA_ITO_GenOSC_v2.mod')
OUTDIR=_os.path.join(ELLIPS_OUT, r'mod_files_v2')
CH=np.load(_os.path.join(ELLIPS_OUT, r'ito_izo_all5.npz'))
GEO={'1':(51.0,3.0,0.00),'2':(51.0,3.0,-0.20),'3':(42.0,2.0,0.18),
     '4':(42.0,2.0,-0.14),'5':(42.0,2.0,0.00)}
NAME={'1':'VendorA_ITO_1','2':'VendorA_ITO_2','3':'VendorB_IZO_1','4':'VendorB_IZO_2','5':'VendorA_ITO_2pctO2'}
GD=0.10   # Drude broadening, eV (fixed)

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
def eps_of(p,E):
    e2=np.maximum(tl_e2(p[2],p[3],p[4],p[5]),0.0)
    e1=p[0]+kk(e2)-p[1]/(EG**2+GD**2)
    e2=e2+p[1]*GD/(EG*(EG**2+GD**2))
    return np.interp(E,EG,e1), np.interp(E,EG,np.maximum(e2,0.0))
def nk_of(p,E):
    e1,e2=eps_of(p,E)
    N=np.sqrt((e1+1j*e2).astype(complex)); N=np.where(N.imag<0,np.conj(N),N)
    return N.real, np.maximum(N.imag,0)

LO=np.array([1.0,0.0,  5.0,3.0,0.05,2.50])
HI=np.array([6.0,4.0,400.0,8.0,8.00,4.50])

def fit_sample(s):
    w,n,k,r=CH[s+'_w'],CH[s+'_n'],CH[s+'_k'],CH[s+'_r']
    m=(w>=360)&(w<=1080)&(r<0.05)
    E=1239.84193/w[m]; e1c=n[m]**2-k[m]**2; e2c=2*n[m]*k[m]; wq=1.0/(0.02+r[m])
    def res(p):
        e1,e2=eps_of(p,E)
        return np.concatenate([wq*(e1-e1c), 3.0*wq*(e2-e2c)])
    best=None
    for Ad0 in [0.05,0.3,0.9,1.8]:
        for Eg0 in [3.0,3.4,3.8]:
            for C0 in [0.5,1.5,3.0]:
                p0=np.clip(np.array([3.0,Ad0,100.,4.8,C0,Eg0]),LO,HI)
                rr=least_squares(res,p0,bounds=(LO,HI),max_nfev=1200)
                if best is None or rr.cost<best.cost: best=rr
    p=best.x
    nf,kf=nk_of(p,E)
    return p, float(np.mean(np.abs(nf-n[m]))), float(np.mean(np.abs(kf-k[m])))

def drude_to_ce(A,Br):
    tau=HBAR_EVFS/Br
    rho=HBAR_EVS**2/(A*EPS0*(tau*1e-15))*100.0
    return rho, tau

def setf(line, idx, val):
    m=re.match(r'^(\t*)(.*)$', line, re.S)
    ind, rest = m.group(1), m.group(2)
    f=rest.split('\t')
    f[idx]=val
    return ind+'\t'.join(f)

def main():
    os.makedirs(OUTDIR,exist_ok=True)
    tpl=open(TEMPLATE,'rb').read().decode('utf-8-sig')
    m=re.search(r"(start_Gen-Osc Fit Parms\r\n)(.*?)(\r\n\t\tend_Gen-Osc Fit Parms)",tpl,re.S)
    head,body,tail=m.group(1),m.group(2),m.group(3)
    tlines=body.split('\r\n')
    assert len(tlines)==11, len(tlines)
    print('template block lines:',len(tlines),'(구조 보존)')
    summary=[]
    for s in '12345':
        p,dn,dk=fit_sample(s)
        d_nm,rg_nm,dth=GEO[s]
        rho,tau=drude_to_ce(p[1],GD)
        L=list(tlines)
        # line0: count stays '2'; field1 = Einf value, field3/4 = bounds
        L[0]=setf(L[0],1,repr(float(p[0]))); L[0]=setf(L[0],3,'0.0'); L[0]=setf(L[0],4,'10.0')
        # line2: Resistivity
        L[2]=setf(L[2],0,repr(float(rho))); L[2]=setf(L[2],2,'1.0E-5'); L[2]=setf(L[2],3,'1.0E-1')
        # line3: Scat. Time  (freeze: not determined by 280-1080nm)
        L[3]=setf(L[3],0,repr(float(tau))); L[3]=setf(L[3],1,'F')
        L[3]=setf(L[3],2,'1.0'); L[3]=setf(L[3],3,'20.0')
        # line5..8: TL Amp, Br, Eo, Eg
        for idx,(v,lo,hi) in zip((5,6,7,8),
                                 ((p[2],'1.0','400.0'),(p[4],'0.05','8.0'),
                                  (p[3],'3.0','8.0'),  (p[5],'2.5','4.5'))):
            L[idx]=setf(L[idx],0,repr(float(v))); L[idx]=setf(L[idx],2,lo); L[idx]=setf(L[idx],3,hi)
        newblock='\r\n'.join(L)
        txt=tpl[:m.start()]+head+newblock+tail+tpl[m.end():]
        # geometry
        txt=re.sub(r"^\t[^\t\n]+(\t[TF]\t-5\.0\t5\.0\tF\t'Angle Offset')",
                   lambda z: '\t'+repr(float(dth))+z.group(1), txt, count=1, flags=re.M)
        txt=re.sub(r"^\t[^\t\n]+(\t[TF]\t0\.0\t500\.0\tF\t'Roughness')",
                   lambda z: '\t'+repr(float(rg_nm*10))+z.group(1), txt, count=1, flags=re.M)
        txt=re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # 2')",
                   lambda z: '\t\t'+repr(float(d_nm*10))+'\tT\t100.0\t2000.0'+z.group(1),
                   txt, count=1, flags=re.M)
        # UV/IR poles OFF: this model's Einf already absorbs the far-UV part
        txt = txt.replace('267.55497728270007\tT\t', '0.0\tF\t', 1)
        txt = txt.replace('11.767447725672241\tT\t', '11.767\tF\t', 1)
        out=os.path.join(OUTDIR,NAME[s]+'_GenOsc.mod')
        open(out,'w',encoding='utf-8-sig',newline='').write(txt)
        nl=txt.count('\r\n')
        summary.append((s,p,rho,tau,dn,dk,nl))
        print('%-16s Einf=%.4f rho=%.4g tau=%.3f | TL(A=%.1f Br=%.3f Eo=%.3f Eg=%.3f) | dn=%.4f dk=%.4f | lines=%d'%(
            NAME[s],p[0],rho,tau,p[2],p[4],p[3],p[5],dn,dk,nl))
    # nk csv from these same params
    wl=np.arange(300.,1101.,2.); E=1239.84193/wl
    cols={s:nk_of(fit_sample(s)[0] if False else summary[i][1],E) for i,s in enumerate('12345')}
    notes=[['# n,k from the SAME parameters as the v2 .mod files (Einf + Drude + Tauc-Lorentz)'],
           ['# Gaussian dropped so the file structure matches CompleteEASE exactly.'],
           ['# Geometry: ITO 51nm/rough 3nm, IZO 42nm/rough 2nm, native oxide 3nm, Si substrate.'],
           ['# Uncertainty on n: +-0.05.  Reliable 360-1080nm.  #5 (2%O2) = do not use.'],['#']]
    with open(os.path.join(OUTDIR,'ITO_IZO_nk_v2.csv'),'w',newline='') as f:
        wtr=csv.writer(f)
        for z in notes: wtr.writerow(z)
        h=['wl_nm']
        for s in '12345': h+=['n_'+NAME[s],'k_'+NAME[s]]
        wtr.writerow(h)
        for i in range(len(wl)):
            row=['%.1f'%wl[i]]
            for s in '12345': row+=['%.4f'%cols[s][0][i],'%.4f'%cols[s][1][i]]
            wtr.writerow(row)
    print('saved ->',OUTDIR)

if __name__=='__main__':
    main()
