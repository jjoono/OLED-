"""T/R/A of MgAg 8 nm film (measured nk #4,#5) on glass, relative to bare glass.
Air / MgAg(8nm) / glass, normal incidence. Compare with pure Ag (#2) 8nm.
"""
import numpy as np, openpyxl, os, csv

XL = r"C:\Users\Junho\Dropbox\Linkstation\Co-work\이선정\Ag15nm, MgAg 8nm_260703.xlsx"
BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
wb = openpyxl.load_workbook(XL, data_only=True); ws = wb["260703"]
# header row 4: C=nm(3) D=n1 E=n2 F=n3 G=n4 H=n5 I=n6 J=k1 K=k2 L=k3 M=k4 N=k5 O=k6
lam=[]; cols={k:[] for k in ["n2","n4","n5","k2","k4","k5"]}
for r in ws.iter_rows(min_row=5, values_only=True):
    if r[2] is None: continue
    lam.append(float(r[2]))
    cols["n2"].append(float(r[4])); cols["n4"].append(float(r[6])); cols["n5"].append(float(r[7]))
    cols["k2"].append(float(r[10])); cols["k4"].append(float(r[12])); cols["k5"].append(float(r[13]))
lam=np.array(lam)
def nk(a,b): return np.array(cols[a])+1j*np.array(cols[b])
MgAg = 0.5*(nk("n4","k4")+nk("n5","k5"))   # average of two MgAg samples
Ag   = nk("n2","k2")

n_air,n_glass=1.0,1.52
def tmm(nl,dl,l):
    M=np.eye(2,dtype=complex)
    for j in range(len(nl)-1):
        a,b=nl[j],nl[j+1]; r=(a-b)/(a+b); t=2*a/(a+b)
        I=np.array([[1,r],[r,1]],dtype=complex)/t
        if j+1<len(nl)-1:
            d=2*np.pi*nl[j+1]*dl[j+1]/l
            M=M@I@np.array([[np.exp(-1j*d),0],[0,np.exp(1j*d)]])
        else: M=M@I
    t_tot=1/M[0,0]; r_tot=M[1,0]/M[0,0]
    T=(np.real(nl[-1])/np.real(nl[0]))*abs(t_tot)**2; R=abs(r_tot)**2
    return T,R

probes=[400,450,500,550,600,650,700]
def at(nkarr,l,d=8.0):
    i=np.argmin(abs(lam-l)); nn=nkarr[i]
    T,R=tmm([n_air,nn,n_glass],[0,d,0],l)
    Tg,_=tmm([n_air,n_glass],[0,0],l)
    return T,R,1-T-R,T/Tg,nn

print("=== MgAg 8 nm (avg of #4,#5) ===")
print(f"{'nm':>4} {'n':>6} {'k':>6} {'T%':>7} {'R%':>7} {'A%':>7} {'T_rel%':>8}")
rows=[]
for l in probes:
    T,R,A,Trel,nn=at(MgAg,l)
    print(f"{l:>4} {nn.real:6.3f} {nn.imag:6.3f} {100*T:7.1f} {100*R:7.1f} {100*A:7.1f} {100*Trel:8.1f}")
    rows.append((l,nn.real,nn.imag,T,R,A,Trel))
print("\n=== pure Ag 8 nm (#2) for comparison ===")
print(f"{'nm':>4} {'n':>6} {'k':>6} {'T%':>7} {'R%':>7} {'A%':>7} {'T_rel%':>8}")
for l in probes:
    T,R,A,Trel,nn=at(Ag,l)
    print(f"{l:>4} {nn.real:6.3f} {nn.imag:6.3f} {100*T:7.1f} {100*R:7.1f} {100*A:7.1f} {100*Trel:8.1f}")

with open(os.path.join(BASE,"MgAg8_TRA.csv"),"w",newline="",encoding="utf-8-sig") as f:
    w=csv.writer(f); w.writerow(["nm","n_MgAg","k_MgAg","T","R","A","T_rel_vs_glass"])
    for l,n,k,T,R,A,Trel in rows: w.writerow([l,f"{n:.4f}",f"{k:.4f}",f"{T:.4f}",f"{R:.4f}",f"{A:.4f}",f"{Trel:.4f}"])
print("\nsaved MgAg8_TRA.csv")
