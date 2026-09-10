"""Cross-validate the ellipsometry result against the absolute T/R campaign.

Their stack: soda-lime glass (1 mm) / HATCN 5 nm / Ag.  Glass is thick, so the
front stack is coherent and the back surface is added incoherently.
Take MY SE-derived Ag n,k and thickness (fitted on Si) and predict their measured
absolute T and R.  If both are right the curves land on each other; the seed
index already agrees (their library l_HATCN n550 = 1.849, my SE fit 1.867).
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import numpy as np, json, csv
import ce_osc as osc, ellipsometry_fit as ef, ag_final as AG, ag_load as L

TR = _os.path.join(ELLIPS_OUT, r'TR260820')
V3 = json.load(open('ag_v3_result.json'))
SEEDP = json.load(open('ag_seed_result.json'))
SP = {k: [r for r in SEEDP[k] if abs(r['d_ox']-2.0) < 1e-9][0]['p'][1:] for k in ('HATCN','MoOx')}

def load_csv(p):
    a = np.genfromtxt(p, delimiter=',', names=True)
    return a

def seedN(p, wl):
    E = 1239.841984/wl; e1, e2 = osc.tauc_lorentz(E, p[1], p[2], p[3], p[4])
    N = np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag < 0, -N, N)

def n_glass(wl):                      # soda-lime, Sellmeier-ish; A~0 in 450-700
    return np.full(len(wl), 1.5230) - 0.0000  # flat is fine over 380-800

def TR_stack(wl, N_ag, d_ag, N_seed, d_seed, theta=0.0):
    """coherent air/Ag/seed/glass, then incoherent glass back surface"""
    amb = np.ones(len(wl)); ng = n_glass(wl).astype(complex)
    # front: air -> Ag -> seed -> glass (semi-infinite)
    rp, rs = ef._tmm(wl, [amb, N_ag, N_seed, ng], [d_ag, d_seed], theta)
    Rf = (np.abs(rp)**2 + np.abs(rs)**2)/2
    # same stack seen from the glass side (for the internal reflection)
    rp2, rs2 = ef._tmm(wl, [ng, N_seed, N_ag, amb], [d_seed, d_ag], theta)
    Rb_in = (np.abs(rp2)**2 + np.abs(rs2)**2)/2
    Tf = 1.0 - Rf - absorbed(wl, N_ag, d_ag, N_seed, d_seed, theta)
    Tf = np.clip(Tf, 0, 1)
    Rg = ((ng.real-1)/(ng.real+1))**2            # glass/air back surface
    Tg = 1.0                                      # transparent window
    denom = 1 - Rb_in*Rg*Tg**2
    T = Tf*Tg*(1-Rg)/denom
    R = Rf + Tf*Tg**2*Rg*(1-Rb_in)/denom * (Tf/np.maximum(1-Rf,1e-9))
    return T, R

def absorbed(wl, N_ag, d_ag, N_seed, d_seed, theta):
    """A of the coherent front stack = 1 - Rf - Tf(into glass); get Tf by energy"""
    amb = np.ones(len(wl)); ng = n_glass(wl).astype(complex)
    rp, rs = ef._tmm(wl, [amb, N_ag, N_seed, ng], [d_ag, d_seed], theta)
    # transmitted power into the glass via the transfer matrix amplitude
    tp, ts = _t_amp(wl, [amb, N_ag, N_seed, ng], [d_ag, d_seed], theta)
    fac = (ng.real/1.0)
    Tf = (np.abs(tp)**2 + np.abs(ts)**2)/2 * fac
    return np.clip(1 - (np.abs(rp)**2+np.abs(rs)**2)/2 - Tf, 0, 1)

def _t_amp(wl, n_layers, d_list, theta_deg):
    """same convention as ef._tmm but returning the transmission amplitudes"""
    th = np.deg2rad(theta_deg); k0 = 2*np.pi/wl
    s0 = n_layers[0]*np.sin(th); N = len(wl)
    q = []
    for ni in n_layers:
        c = np.sqrt((ni**2 - s0**2).astype(complex)); q.append(np.where(c.imag < 0, -c, c))
    def eye(): 
        M = np.zeros((N,2,2), complex); M[:,0,0]=1; M[:,1,1]=1; return M
    mm = lambda A,B: np.einsum('nij,njk->nik', A, B)
    def I_(r,t):
        M=np.zeros((N,2,2),complex); M[:,0,0]=1/t; M[:,0,1]=r/t; M[:,1,0]=r/t; M[:,1,1]=1/t; return M
    def Ip(ni,nj,qi,qj):
        d=nj**2*qi+ni**2*qj; return I_((nj**2*qi-ni**2*qj)/d, 2*ni*nj*qi/d)
    def Is(qi,qj):
        d=qi+qj; return I_((qi-qj)/d, 2*qi/d)
    def Pr(qj,d):
        b=k0*qj*d; M=np.zeros((N,2,2),complex)
        M[:,0,0]=np.exp(-1j*b); M[:,1,1]=np.exp(1j*b); return M
    Mp,Ms = eye(),eye()
    for i,d in enumerate(d_list):
        ni,nj,qi,qj = n_layers[i],n_layers[i+1],q[i],q[i+1]
        Mp = mm(mm(Mp, Ip(ni,nj,qi,qj)), Pr(qj,d)); Ms = mm(mm(Ms, Is(qi,qj)), Pr(qj,d))
    ni,nj,qi,qj = n_layers[-2],n_layers[-1],q[-2],q[-1]
    Mp = mm(Mp, Ip(ni,nj,qi,qj)); Ms = mm(Ms, Is(qi,qj))
    return np.conj(1/Mp[:,0,0]), np.conj(1/Ms[:,0,0])

D = load_csv(TR+r'\ALL_SAMPLES_TRA.csv')
wl_m = D['wavelength_nm']
sel = (wl_m>=400)&(wl_m<=800); wl = np.sort(wl_m[sel])
MAP = {4:'1-6',5:'1-7',6:'1-8',7:'2-5',8:'2-6',10:'2-7',12:'2-8'}
print('HATCN series: my SE model vs their measured absolute T / R')
print('%-6s %-7s | %-19s | %-19s'%('Ag nm','lambda','T meas / T model','R meas / R model'))
for nom in (5,7,8,12):
    sh = MAP[nom]; p=np.array(V3[sh]['p'])
    Ns = seedN(SP['HATCN'], wl); Na = AG.agN(p, wl)
    T,R = TR_stack(wl, Na, p[0], Ns, 5.0)
    Tm = np.interp(wl, np.sort(wl_m[sel]), D['HATCN5_Ag%d_T_pct'%nom][sel][np.argsort(wl_m[sel])])/100
    Rm = np.interp(wl, np.sort(wl_m[sel]), D['HATCN5_Ag%d_Rcorr_pct'%nom][sel][np.argsort(wl_m[sel])])/100
    for t in (450,550,650,750):
        i=np.argmin(abs(wl-t))
        print('%-6d %-7.0f | %6.1f%% / %6.1f%%    | %6.1f%% / %6.1f%%'%(
            nom,wl[i],100*Tm[i],100*T[i],100*Rm[i],100*R[i]))
    m=(wl>=430)&(wl<=780)
    print('   -> mean |dT| = %.2f %%p ,  mean |dR| = %.2f %%p  (d_Ag = %.2f nm from SE)'%(
        100*np.mean(np.abs(T[m]-Tm[m])), 100*np.mean(np.abs(R[m]-Rm[m])), p[0]))
