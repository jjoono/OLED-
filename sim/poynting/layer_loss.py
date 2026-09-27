"""Layer-resolved absorption from the net Poynting flux at every interface (transfer matrix, s and p),
for light incident from the substrate on ITO / organics 420 nm (n 1.8) / reflector 100 nm / air.
Weighting: Ir(ppy)2acac spectrum (photon number, 400-700 nm) x Lambertian (cos.sin) over substrate angles.
Absorption of layer j = S_z(entering j) - S_z(leaving j); checked against 1 - R - T."""
import numpy as np, sys, csv
sys.path.insert(0, '/home/user/OLED-/sim/design_rule4'); import fig2d as F
LAM = np.arange(400.0, 701.0, 5.0); sel = np.isin(F.LAM, LAM)
TH = np.radians(np.arange(90) + 0.5); g = np.cos(TH)*np.sin(TH); g /= g.sum()
S = np.clip(np.interp(LAM, F.LAM, F.GREEN), 0, None)*LAM; S /= S.sum()
def absorb(ns, layers, nexit):
    """returns per-layer absorbed fraction (Nl x Nth x Nlayers), R, T (unpolarised)"""
    nl = [np.full(LAM.shape, ns, complex)] + [m for m, _ in layers] + [nexit]
    dl = [0.0] + [d for _, d in layers] + [0.0]
    k0 = 2*np.pi/LAM[:, None]; kx = ns*np.sin(TH)[None, :]*k0
    kz = [np.sqrt((n[:, None]*k0)**2 - kx**2 + 0j) for n in nl]
    kz = [np.where(np.imag(q) < 0, -q, q) for q in kz]
    out = {}
    for pol in 'sp':
        eta = kz if pol == 's' else [kz[j]/nl[j][:, None]**2 for j in range(len(nl))]
        N = len(nl)
        # amplitudes (forward a, backward b) at the left side of each layer, backward propagation from the exit
        a = [None]*N; b = [None]*N
        a[-1] = np.ones_like(kx, complex); b[-1] = np.zeros_like(kx, complex)
        for j in range(N-2, -1, -1):
            # field at right side of layer j equals field at left of j+1
            aj1, bj1 = a[j+1], b[j+1]
            E = aj1 + bj1; H = eta[j+1]*(aj1 - bj1)
            ar = 0.5*(E + H/eta[j]); br = 0.5*(E - H/eta[j])
            ph = np.exp(1j*kz[j]*dl[j])
            a[j] = ar/ph; b[j] = br*ph              # to the left side of layer j
        inc = np.abs(a[0])**2*np.real(eta[0])
        def Sz(j, side):                            # net flux at left (0) / right (1) side of layer j
            ph = np.exp(1j*kz[j]*dl[j]) if side else 1.0
            aa, bb = a[j]*ph, b[j]/ph
            if pol == 's': return np.real((aa+bb)*np.conj(eta[j]*(aa-bb)))
            return np.real(np.conj(aa+bb)*eta[j]*(aa-bb))
        R = np.abs(b[0]/a[0])**2
        T = Sz(N-1, 0)/inc
        A = np.stack([(Sz(j, 0) - Sz(j, 1))/inc for j in range(1, N-1)], -1)
        out[pol] = (A, R, T)
    return tuple(0.5*(out['s'][i] + out['p'][i]) for i in range(3))
if __name__ == '__main__':
    ITO, ORG, AG, AL = F.ITO[sel], F.ORG[sel], F.AG[sel], F.AL[sel]; AIR = np.ones(LAM.shape, complex)
    rows = [['n_sub', 'reflector', 'd_ITO_nm', 'A_ITO', 'A_organics', 'A_reflector', "A' = 1-R", 'T', 'closure 1-R-T-sumA']]
    for ns in (1.5, 1.8):
        for name, M in (('Al', AL), ('Ag', AG)):
            for d in range(30, 201, 10):
                A, R, T = absorb(ns, [(ITO, float(d)), (ORG, 420.0), (M, 100.0)], AIR)
                w = lambda X: float(S @ (X @ g))
                a = [w(A[..., j]) for j in range(3)]; r = w(R); t = w(T)
                rows.append([ns, name, d] + a + [1 - r, t, 1 - r - t - sum(a)])
    with open('/home/user/OLED-/sim/poynting/layer_loss.csv', 'w', newline='') as f: csv.writer(f).writerows(rows)
    for r in rows:
        if r[0] == 'n_sub' or r[2] in (50, 150): print(r if r[0]=='n_sub' else [r[0], r[1], r[2]] + [round(x, 4) for x in r[3:]])
