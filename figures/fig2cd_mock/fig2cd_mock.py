"""
Fig. 2(c)(d) draft: round-trip loss of the OLED stack seen by recycled light.

(c) Spectrally weighted round-trip loss A'(theta) = 1 - R_LED(theta) for three
    reflector systems, plane wave incident from the substrate (n = 1.77).
(d) Layer-resolved net Poynting flux S_z(z) for plane-wave incidence, averaged
    over the returned-light angular distribution and the EL spectrum. The value
    in the substrate equals A'; each step is the loss in that layer; a non-zero
    tail in the rear medium is leakage.

Illustrative only: textbook n,k, simplified stacks. Replace with measured n,k,
the real stack, and the MLA-returned distribution P_sub * BSDF_R.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---------------------------------------------------------------- optical constants (textbook-ish)
WL = np.arange(440, 641, 5.0)            # nm, covers the green EL band
def interp(tab):
    l, n, k = np.array(tab).T
    return np.interp(WL, l, n) + 1j*np.interp(WL, l, k)

Ag = interp([(400,0.05,2.07),(450,0.04,2.66),(500,0.05,3.13),(550,0.055,3.62),
             (600,0.06,4.15),(650,0.065,4.5),(700,0.07,4.84)])
Al = interp([(400,0.49,4.86),(450,0.62,5.47),(500,0.77,6.08),(550,0.96,6.69),
             (600,1.20,7.26),(650,1.47,7.79),(700,1.83,8.31)])
K_TCO = 0.02                                   # as-deposited IZO band (Fig. 3c abscissa)
IZO = 2.0 + 0.0j*WL + 1j*K_TCO
ORG = 1.80 + 0j*WL
ZnS = 2.35 + 0j*WL
LiF = 1.39 + 0j*WL
AIR = 1.0 + 0j*WL
GLASS = 1.50 + 0j*WL
N_SUB = 1.77                                    # incidence medium (substrate / MLA)

lam0 = 540.0
d_H, d_L = lam0/(4*2.35), lam0/(4*1.39)

def stack(kind, rear="air"):
    """list of (name, n(WL) complex array or scalar, thickness nm); first/last semi-infinite"""
    common = [("substrate", N_SUB+0j*WL, None), ("TCO", IZO, 50.0), ("organics", ORG, 300.0)]
    if kind == "Al":
        return common + [("Al", Al, 100.0), ("air", AIR, None)]
    if kind == "Ag":
        return common + [("Ag", Ag, 100.0), ("air", AIR, None)]
    if kind == "DBR":
        dbr = []
        for p in range(4):
            dbr += [("ZnS", ZnS, d_H), ("LiF", LiF, d_L)]
        dbr += [("ZnS", ZnS, d_H)]
        rear_n = AIR if rear == "air" else GLASS
        return common + [("TCO2", IZO, 50.0)] + dbr + [(rear, rear_n, None)]
    raise ValueError

# ---------------------------------------------------------------- TMM with layer-resolved S_z
def solve(layers, theta_deg, pol):
    """Returns R, T, and a function giving normalized S_z on a depth grid.
    Amplitudes A_j (forward), B_j (backward) referenced to each layer's front boundary."""
    th = np.deg2rad(theta_deg)
    k0 = 2*np.pi/WL                                   # 1/nm, array over WL
    n_list = [L[1] for L in layers]
    d_list = [L[2] for L in layers]
    kx = k0*N_SUB*np.sin(th)
    eps = [n**2 for n in n_list]
    q = [np.sqrt(e*k0**2 - kx**2 + 0j) for e in eps]  # principal branch: Im>=0
    # "impedance" quantity: s -> q ; p -> q/eps
    p = [qj if pol == "s" else qj/e for qj, e in zip(q, eps)]
    N = len(layers)
    # backward recursion of effective reflection rho_j = B_j/A_j at front of layer j
    rho = [None]*N
    rho[N-1] = np.zeros_like(WL, dtype=complex)
    r = [None]*(N-1)
    t = [None]*(N-1)
    for j in range(N-1):
        r[j] = (p[j]-p[j+1])/(p[j]+p[j+1])
        t[j] = 2*p[j]/(p[j]+p[j+1])
    for j in range(N-2, -1, -1):
        ph = np.exp(2j*q[j+1]*d_list[j+1]) if d_list[j+1] is not None else 1.0
        rho_back = rho[j+1]*ph  # B/A referenced to the BACK of layer j+1? -> careful:
        # rho[j+1] is B_{j+1}/A_{j+1} at the FRONT of layer j+1. The Airy formula at the
        # interface j|j+1 needs the ratio (b/a) just after the interface = rho[j+1]. Then
        # rho[j] = (r + rho[j+1])/(1 + r rho[j+1]) * exp(2 i q_j d_j).
        ratio = (r[j] + rho[j+1])/(1 + r[j]*rho[j+1])
        rho[j] = ratio*(np.exp(2j*q[j]*d_list[j]) if d_list[j] is not None else 1.0)
    # forward pass
    A = [None]*N; B = [None]*N
    A[0] = np.ones_like(WL, dtype=complex); B[0] = rho[0]
    for j in range(N-1):
        a = A[j]*(np.exp(1j*q[j]*d_list[j]) if d_list[j] is not None else 1.0)
        A[j+1] = t[j]*a/(1 + r[j]*rho[j+1])
        B[j+1] = rho[j+1]*A[j+1]
    R = np.abs(B[0])**2
    T = np.real(p[N-1])/np.real(p[0])*np.abs(A[N-1])**2
    def Sz_profile(zgrid_per_layer):
        """zgrid_per_layer: list of arrays of local depth (0..d) per layer. Returns S_z/S_inc."""
        out = []
        S_inc = np.real(p[0])  # |A0|^2 = 1
        for j, zz in enumerate(zgrid_per_layer):
            f = A[j][:, None]*np.exp(1j*q[j][:, None]*zz[None, :])
            b = B[j][:, None]*np.exp(-1j*q[j][:, None]*zz[None, :])
            S = np.real(p[j][:, None]*(f - b)*np.conj(f + b)) if pol == "p" else \
                np.real(np.conj(p[j][:, None]*(f - b))*(f + b))
            # both forms equal Re[p (f-b) conj(f+b)] for s; for p the same expression
            # with p = q/eps; sign conventions identical -> use one:
            S = np.real(p[j][:, None]*(f - b)*np.conj(f + b))
            out.append(S/S_inc[:, None])
        return out
    return R, T, Sz_profile

# ---------------------------------------------------------------- weights
EL = np.exp(-0.5*((WL-530.0)/25.5)**2); EL /= EL.sum()          # FWHM ~ 60 nm
thetas = np.arange(0, 89.5, 1.0)
w_th = np.sin(np.deg2rad(thetas))*np.cos(np.deg2rad(thetas)); w_th /= w_th.sum()  # Lambertian in substrate (proxy for MLA return)

def A_prime_vs_theta(layers):
    out = []
    for th in thetas:
        Rs, _, _ = solve(layers, th, "s"); Rp, _, _ = solve(layers, th, "p")
        out.append(np.sum(EL*(1 - 0.5*(Rs+Rp))))
    return np.array(out)

def weighted_profile(layers, pad=120.0, dz=1.0):
    d_list = [L[2] for L in layers]
    grids = []
    for j, d in enumerate(d_list):
        if d is None:
            grids.append(np.arange(-pad, 0, dz) if j == 0 else np.arange(0, pad, dz))
        else:
            grids.append(np.arange(0, d+1e-9, dz))
    acc = [np.zeros(len(g)) for g in grids]
    for th, w in zip(thetas, w_th):
        for pol in ("s", "p"):
            _, _, prof = solve(layers, th, pol)
            S = prof(grids)
            for j in range(len(grids)):
                acc[j] += 0.5*w*np.sum(EL[:, None]*S[j], axis=0)
    # absolute depth axis
    z0 = 0.0; zabs = []
    for j, g in enumerate(grids):
        if j == 0: zabs.append(g)            # substrate: negative local depth (already -pad..0)
        else:
            zabs.append(z0 + g); z0 += (d_list[j] if d_list[j] is not None else 0.0)
    return zabs, acc, d_list

# ---------------------------------------------------------------- compute
systems = [("Al / TCO",  stack("Al"),  "#2a78d6"),
           ("Ag / TCO",  stack("Ag"),  "#eb6834"),
           ("DBR+TCO / TCO (air-backed)", stack("DBR","air"), "#1baf7a")]
dbr_glass = stack("DBR","glass")

Ap = {name: A_prime_vs_theta(L) for name, L, _ in systems}
Ap_glass = A_prime_vs_theta(dbr_glass)
Ap_w = {name: np.sum(w_th*Ap[name]) for name in Ap}

profiles = {name: weighted_profile(L) for name, L, _ in systems}

# per-layer loss table
def layer_losses(name):
    zabs, acc, d_list = profiles[name]
    rows = []
    layers = dict(systems_map)[name]
    names = [L[0] for L in layers]
    for j in range(1, len(acc)-1):
        rows.append((names[j], acc[j][0] - acc[j][-1]))
    return acc[0][0], rows, acc[-1][0]
systems_map = [(n, L) for n, L, _ in systems]

# ---------------------------------------------------------------- figure
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "font.family": "DejaVu Sans"})
fig = plt.figure(figsize=(11.0, 6.4), facecolor="white")
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.55, wspace=0.35,
                      left=0.07, right=0.98, top=0.93, bottom=0.10)
axc = fig.add_subplot(gs[0, :])
axd = [fig.add_subplot(gs[1, i]) for i in range(3)]

# ---- (c)
axc.fill_between(thetas, 0, w_th/w_th.max()*0.12, color="#e6e6e3", lw=0, label="returned-light weight (a.u.)")
for name, L, col in systems:
    axc.plot(thetas, Ap[name], color=col, lw=2.0, label=name)
axc.plot(thetas, Ap_glass, color="#1baf7a", lw=1.6, ls="--", label="DBR+TCO / TCO (glass-backed)")
tir_air = np.degrees(np.arcsin(1/N_SUB)); tir_gl = np.degrees(np.arcsin(1.5/N_SUB))
for x, lab in [(tir_air, "TIR, air-backed"), (tir_gl, "TIR, glass-backed")]:
    axc.axvline(x, color="#9a9a96", lw=0.8, ls=":")
    axc.text(x+0.6, 0.30, lab, rotation=90, va="top", ha="left", fontsize=7.5, color="#52514e")
summary = ("Angle-weighted round-trip loss ⟨A′⟩\n"
           f"Al / TCO: {Ap_w['Al / TCO']*100:.1f}%\n"
           f"Ag / TCO: {Ap_w['Ag / TCO']*100:.1f}%\n"
           f"DBR+TCO / TCO, air-backed: {Ap_w['DBR+TCO / TCO (air-backed)']*100:.1f}%\n"
           f"DBR+TCO / TCO, glass-backed: {np.sum(w_th*Ap_glass)*100:.1f}%")
axc.text(0.012, 0.97, summary, transform=axc.transAxes, ha="left", va="top", fontsize=8,
         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#d0d0cc", lw=0.6))
axc.text(46, 0.77, "leaks into glass\n(no TIR below 58°)", fontsize=7.5, color="#52514e", ha="center", va="top")
axc.set_xlim(0, 90); axc.set_ylim(0, 0.8)
axc.set_xlabel("Angle of incidence in substrate, θ (deg)")
axc.set_ylabel("Round-trip loss  A′(θ) = 1 − R_LED")
axc.set_title("(c)  Round-trip loss of the OLED stack seen by recycled light (EL-weighted, s/p average)", loc="left", fontsize=10)
axc.legend(frameon=False, fontsize=8, loc="upper right", ncol=1)
for s in ("top", "right"): axc.spines[s].set_visible(False)
axc.grid(axis="y", color="#e6e6e3", lw=0.6)

# ---- (d)
band_col = {"substrate": "#dfe8f3", "TCO": "#f6dfd0", "TCO2": "#f6dfd0", "organics": "#eef2e6",
            "Al": "#b9b9b6", "Ag": "#d9d9d6", "ZnS": "#e3e3f5", "LiF": "#f4f4fb", "air": "#ffffff", "glass": "#eef4fa"}
for ax, (name, L, col) in zip(axd, systems):
    zabs, acc, d_list = profiles[name]
    names = [l[0] for l in L]
    # layer bands
    for j in range(len(zabs)):
        z0, z1 = zabs[j][0], zabs[j][-1]
        ax.axvspan(z0, z1, color=band_col.get(names[j], "#eeeeee"), lw=0)
    zz = np.concatenate(zabs); ss = np.concatenate(acc)
    ax.plot(zz, ss, color=col, lw=2.0)
    A0 = acc[0][0]; T = acc[-1][0]
    ax.axhline(A0, color="#9a9a96", lw=0.7, ls=":")
    ax.text(zabs[0][0]+3, A0+0.006, f"A′ = {A0*100:.1f}%", fontsize=8, va="bottom")
    # per-layer annotations (merge DBR layers)
    losses = {}
    for j in range(1, len(acc)-1):
        key = "DBR" if names[j] in ("ZnS", "LiF") else names[j]
        losses[key] = losses.get(key, 0.0) + (acc[j][0]-acc[j][-1])
    txt = "\n".join(f"{k}: {v*100:.1f}%" for k, v in losses.items() if abs(v) > 5e-4)
    if T > 5e-4: txt += f"\nleakage: {T*100:.1f}%"
    ax.text(0.98, 0.72, txt, transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d0d0cc", lw=0.6))
    ax.set_ylim(0, max(0.30, A0*1.25)); ax.set_xlim(zz[0], zz[-1])
    ax.set_xlabel("Depth from substrate/TCO interface (nm)")
    ax.set_title(f"(d) {name}", loc="left", fontsize=9.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    # layer labels on top
    for j in range(len(zabs)):
        if names[j] in ("ZnS", "LiF"): continue
        zm = 0.5*(zabs[j][0]+zabs[j][-1])
        ax.text(zm, ax.get_ylim()[1]*0.985, names[j], ha="center", va="top", fontsize=6.5, color="#52514e")
    if "DBR" in name:
        zdbr = [zabs[j] for j in range(len(zabs)) if names[j] in ("ZnS", "LiF")]
        zm = 0.5*(zdbr[0][0]+zdbr[-1][-1])
        ax.text(zm, ax.get_ylim()[1]*0.985, "DBR (ZnS/LiF ×4.5)", ha="center", va="top", fontsize=6.5, color="#52514e")
axd[0].set_ylabel("Net flux S_z / returned flux")

fig.text(0.07, 0.015, "Draft with textbook n,k; substrate n=1.77, TCO 50 nm (k=0.02), lossless organics 300 nm, "
         "reflector 100 nm or 4.5-pair ZnS/LiF; Lambertian return proxy, EL Gaussian 530/60 nm.",
         fontsize=7.5, color="#52514e")
out = "/home/user/OLED-/figures/fig2cd_mock/fig2cd_mock.png"
fig.savefig(out, dpi=170)
print("saved", out)
for name in Ap:
    A0, rows, T = layer_losses(name)
    print(f"{name}: A'={A0*100:.2f}%  T={T*100:.2f}%  ", [(r[0], round(r[1]*100, 2)) for r in rows])
print("A'(0deg):", {k: round(v[0]*100, 2) for k, v in Ap.items()}, " glass-backed DBR <A'>:", round(np.sum(w_th*Ap_glass)*100, 2))
