"""Round-trip check: parse the generated .mod files back, convert
CompleteEASE parameters to the internal convention, rebuild n,k and compare
with the target dispersion (genosc_params.json + delivered CSV)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, json, math, csv, os
import numpy as np

HBAR_EVS=6.582119569e-16; HBAR_EVFS=0.6582119569; EPS0=8.8541878128e-12
OUTDIR=_os.path.join(ELLIPS_OUT, r'mod_files')
NAMEMAP={'1':'VendorA_ITO_1','2':'VendorA_ITO_2','3':'VendorB_IZO_1','4':'VendorB_IZO_2','5':'VendorA_ITO_2pctO2'}

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

def nk_from(einf,Ad,Gd,tlp,gp,E):
    e2=np.maximum(tl_e2(*tlp)+gau_e2(*gp),0.0)
    e1=einf+kk(e2)-Ad/(EG**2+Gd**2)
    e2=e2+Ad*Gd/(EG*(EG**2+Gd**2))
    N=np.sqrt((np.interp(E,EG,e1)+1j*np.interp(E,EG,np.maximum(e2,0))).astype(complex))
    N=np.where(N.imag<0,np.conj(N),N)
    return N.real, np.maximum(N.imag,0)

def parse_mod(path):
    t=open(path,'rb').read().decode('utf-8-sig')
    blk=re.search(r"start_Gen-Osc Fit Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Fit Parms",t,re.S).group(1)
    lines=[l.strip() for l in blk.split('\r\n') if l.strip()]
    einf=float(lines[0].split('\t')[1])
    vals={}; order=[]
    i=1
    while i<len(lines):
        L=lines[i]
        if L.startswith("'"):
            typ=L.strip().strip("'"); order.append(typ); vals[typ]=[]; i+=1
            while i<len(lines) and not lines[i].startswith("'") and lines[i]!='F':
                vals[typ].append(float(lines[i].split('\t')[0])); i+=1
        else: i+=1
    d=float(re.search(r"^\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Thickness # 2'",t,re.M).group(1))
    rough=float(re.search(r"^\t([\d.eE+-]+)\t[TF]\t0\.0\t500\.0\tF\t'Roughness'",t,re.M).group(1))
    dth=float(re.search(r"^\t([-\d.eE+]+)\t[TF]\t-5\.0\t5\.0\tF\t'Angle Offset'",t,re.M).group(1))
    return einf,vals,order,d,rough,dth

# target CSV (delivered)
W=[];COLS=None;DATA=[]
for row in csv.reader(open(_os.path.join(ELLIPS_OUT, r'ITO_IZO_final_nk.csv'))):
    if not row or row[0].startswith('#'): continue
    if COLS is None: COLS=row; continue
    DATA.append([float(x) for x in row])
DATA=np.array(DATA); wl_csv=DATA[:,0]

R=json.load(open(_os.path.join(ELLIPS_OUT, r'genosc_params.json')))
Etest=1239.84193/np.array([400.,450.,500.,550.,633.,700.,800.,1000.])
print('%-16s %9s %9s | %9s %9s | geometry'%('sample','max|dn|*','max|dk|*','max|dn|CSV','max|dk|CSV'))
ok=True
for s in '12345':
    p=parse_mod(os.path.join(OUTDIR,NAMEMAP[s]+'_GenOsc.mod'))
    einf,vals,order,d_A,rough_A,dth=p
    rho,tau=vals['Drude(RT)']
    Ad=HBAR_EVS**2/((rho/100.0)*EPS0*(tau*1e-15)); Gd=HBAR_EVFS/tau
    tA,tBr,tEo,tEg=vals['Tauc-Lorentz']
    gA,gBr,gEn=vals['Gaussian']
    sig=gBr/(2*math.sqrt(math.log(2)))
    n_mod,k_mod=nk_from(einf,Ad,Gd,(tA,tEo,tBr,tEg),(gA,gEn,sig),Etest)
    # reference: direct from stored params
    q=R[s]['p']
    n_ref,k_ref=nk_from(q[0],q[1],0.10,(q[2],q[3],q[4],q[5]),(q[6],q[7],q[8]),Etest)
    dn=np.max(np.abs(n_mod-n_ref)); dk=np.max(np.abs(k_mod-k_ref))
    i=[int(np.argmin(np.abs(wl_csv-w))) for w in [400,450,500,550,633,700,800,1000]]
    ncsv=DATA[i,1+2*int(s)-2]; kcsv=DATA[i,2+2*int(s)-2]
    dnc=np.max(np.abs(n_mod-ncsv)); dkc=np.max(np.abs(k_mod-kcsv))
    flag='' if (dn<1e-6 and dk<1e-6) else '  <-- MISMATCH'
    if dn>1e-6 or dk>1e-6: ok=False
    print('%-16s %9.2e %9.2e | %9.4f %9.4f | d=%.0fA rough=%.0fA dth=%+.2f%s'%(
        NAMEMAP[s],dn,dk,dnc,dkc,d_A,rough_A,dth,flag))
print()
print('* = .mod 역파싱 값 vs 저장된 파라미터 (변환 정확도, 0이어야 정상)')
print('CSV 열은 별도 피팅본이라 소폭 차이는 정상')
print('RESULT:', 'PASS - 변환 무손실' if ok else 'FAIL')
