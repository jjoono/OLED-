"""
Fig. 2(c)(d) draft (main text): round-trip loss seen by recycled light, Al vs Ag,
with the transparent-electrode extinction k as the second variable.
Supplementary figure: DBR behind a TCO, air- vs glass-backed.

Plane-wave incidence from the substrate (n = 1.77); transfer matrix with
layer-resolved net Poynting flux. Textbook n,k, simplified stacks: illustrative.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WL = np.arange(440, 641, 5.0)
def interp(tab):
    l, n, k = np.array(tab).T
    return np.interp(WL, l, n) + 1j*np.interp(WL, l, k)
Ag = interp([(400,0.05,2.07),(450,0.04,2.66),(500,0.05,3.13),(550,0.055,3.62),(600,0.06,4.15),(650,0.065,4.5),(700,0.07,4.84)])
Al = interp([(400,0.49,4.86),(450,0.62,5.47),(500,0.77,6.08),(550,0.96,6.69),(600,1.20,7.26),(650,1.47,7.79),(700,1.83,8.31)])
ORG, ZnS, LiF, AIR, GLASS = (1.80+0j*WL), (2.35+0j*WL), (1.39+0j*WL), (1.0+0j*WL), (1.50+0j*WL)
N_SUB = 1.77
def TCO(k): return 2.0 + 0j*WL + 1j*k
lam0 = 540.0; d_H, d_L = lam0/(4*2.35), lam0/(4*1.39)

def stack(kind, k_tco=0.02, rear="air"):
    common = [("substrate", N_SUB+0j*WL, None), ("TCO", TCO(k_tco), 50.0), ("organics", ORG, 300.0)]
    if kind == "Al": return common + [("Al", Al, 100.0), ("air", AIR, None)]
    if kind == "Ag": return common + [("Ag", Ag, 100.0), ("air", AIR, None)]
    if kind == "DBR":
        dbr = []
        for _ in range(4): dbr += [("ZnS", ZnS, d_H), ("LiF", LiF, d_L)]
        dbr += [("ZnS", ZnS, d_H)]
        return common + [("TCO2", TCO(k_tco), 50.0)] + dbr + [(rear, AIR if rear == "air" else GLASS, None)]

def solve(layers, theta_deg, pol):
    th = np.deg2rad(theta_deg); k0 = 2*np.pi/WL
    n_list = [L[1] for L in layers]; d_list = [L[2] for L in layers]
    kx = k0*N_SUB*np.sin(th); eps = [n**2 for n in n_list]
    q = [np.sqrt(e*k0**2 - kx**2 + 0j) for e in eps]
    p = [qj if pol == "s" else qj/e for qj, e in zip(q, eps)]
    N = len(layers)
    r = [(p[j]-p[j+1])/(p[j]+p[j+1]) for j in range(N-1)]
    t = [2*p[j]/(p[j]+p[j+1]) for j in range(N-1)]
    rho = [None]*N; rho[N-1] = np.zeros_like(WL, dtype=complex)
    for j in range(N-2, -1, -1):
        ratio = (r[j] + rho[j+1])/(1 + r[j]*rho[j+1])
        rho[j] = ratio*(np.exp(2j*q[j]*d_list[j]) if d_list[j] is not None else 1.0)
    A = [None]*N; B = [None]*N
    A[0] = np.ones_like(WL, dtype=complex); B[0] = rho[0]
    for j in range(N-1):
        a = A[j]*(np.exp(1j*q[j]*d_list[j]) if d_list[j] is not None else 1.0)
        A[j+1] = t[j]*a/(1 + r[j]*rho[j+1]); B[j+1] = rho[j+1]*A[j+1]
    R = np.abs(B[0])**2
    T = np.real(p[N-1])/np.real(p[0])*np.abs(A[N-1])**2
    def Sz_profile(grids):
        out = []; S_inc = np.real(p[0])
        for j, zz in enumerate(grids):
            f = A[j][:, None]*np.exp(1j*q[j][:, None]*zz[None, :])
            b = B[j][:, None]*np.exp(-1j*q[j][:, None]*zz[None, :])
            out.append(np.real(p[j][:, None]*(f - b)*np.conj(f + b))/S_inc[:, None])
        return out
    return R, T, Sz_profile

EL = np.exp(-0.5*((WL-530.0)/25.5)**2); EL /= EL.sum()
thetas = np.arange(0, 89.5, 1.0)
w_th = np.sin(np.deg2rad(thetas))*np.cos(np.deg2rad(thetas)); w_th /= w_th.sum()

def A_prime_vs_theta(layers):
    out = []
    for th in thetas:
        Rs, _, _ = solve(layers, th, "s"); Rp, _, _ = solve(layers, th, "p")
        out.append(np.sum(EL*(1 - 0.5*(Rs+Rp))))
    return np.array(out)

def weighted_profile(layers, pad=120.0, dz=1.0):
    d_list = [L[2] for L in layers]; grids = []
    for j, d in enumerate(d_list):
        grids.append((np.arange(-pad, 0, dz) if j == 0 else np.arange(0, pad, dz)) if d is None else np.arange(0, d+1e-9, dz))
    acc = [np.zeros(len(g)) for g in grids]
    for th, w in zip(thetas, w_th):
        for pol in ("s", "p"):
            _, _, prof = solve(layers, th, pol); S = prof(grids)
            for j in range(len(grids)): acc[j] += 0.5*w*np.sum(EL[:, None]*S[j], axis=0)
    z0 = 0.0; zabs = []
    for j, g in enumerate(grids):
        if j == 0: zabs.append(g)
        else: zabs.append(z0 + g); z0 += (d_list[j] or 0.0)
    return zabs, acc

BAND = {"substrate": "#dfe8f3", "TCO": "#f6dfd0", "TCO2": "#f6dfd0", "organics": "#eef2e6", "Al": "#b9b9b6",
        "Ag": "#d9d9d6", "ZnS": "#e3e3f5", "LiF": "#f4f4fb", "air": "#ffffff", "glass": "#eef4fa"}

def draw_profile(ax, layers, col, title, ymax=0.30):
    zabs, acc = weighted_profile(layers); names = [l[0] for l in layers]
    for j in range(len(zabs)): ax.axvspan(zabs[j][0], zabs[j][-1], color=BAND.get(names[j], "#eee"), lw=0)
    zz = np.concatenate(zabs); ss = np.concatenate(acc); ax.plot(zz, ss, color=col, lw=2.0)
    A0 = acc[0][0]; T = acc[-1][0]
    ax.axhline(A0, color="#9a9a96", lw=0.7, ls=":"); ax.text(zabs[0][0]+3, A0+0.006, f"A′ = {A0*100:.1f}%", fontsize=8, va="bottom")
    losses = {}
    for j in range(1, len(acc)-1):
        key = "DBR" if names[j] in ("ZnS", "LiF") else names[j]
        losses[key] = losses.get(key, 0.0) + (acc[j][0]-acc[j][-1])
    txt = "\n".join(f"{k}: {v*100:.1f}%" for k, v in losses.items() if abs(v) > 5e-4)
    if T > 5e-4: txt += f"\nleakage: {T*100:.1f}%"
    ax.text(0.98, 0.72, txt, transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#d0d0cc", lw=0.6))
    ax.set_ylim(0, ymax); ax.set_xlim(zz[0], zz[-1]); ax.set_title(title, loc="left", fontsize=9.5)
    ax.set_xlabel("Depth from substrate/TCO interface (nm)")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    for j in range(len(zabs)):
        if names[j] in ("ZnS", "LiF"): continue
        ax.text(0.5*(zabs[j][0]+zabs[j][-1]), ymax*0.985, names[j], ha="center", va="top", fontsize=6.5, color="#52514e")
    if any(n in ("ZnS", "LiF") for n in names):
        zd = [zabs[j] for j in range(len(zabs)) if names[j] in ("ZnS", "LiF")]
        ax.text(0.5*(zd[0][0]+zd[-1][-1]), ymax*0.985, "DBR (ZnS/LiF ×4.5)", ha="center", va="top", fontsize=6.5, color="#52514e")
    return A0

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.6, "font.family": "DejaVu Sans"})
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
K_HI, K_LO = 0.02, 0.005   # as-deposited IZO vs annealed ITO (Fig. 3c bands)

# ============================================================ MAIN FIGURE: Al vs Ag
fig = plt.figure(figsize=(11.0, 6.4), facecolor="white")
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.35, left=0.07, right=0.98, top=0.93, bottom=0.10)
axc = fig.add_subplot(gs[0, :]); axd = [fig.add_subplot(gs[1, i]) for i in range(3)]

curves = [("Al reflector, TCO k = 0.02", stack("Al", K_HI), BLUE, "-"),
          ("Al reflector, TCO k = 0.005", stack("Al", K_LO), BLUE, "--"),
          ("Ag reflector, TCO k = 0.02", stack("Ag", K_HI), ORANGE, "-"),
          ("Ag reflector, TCO k = 0.005", stack("Ag", K_LO), ORANGE, "--")]
axc.fill_between(thetas, 0, w_th/w_th.max()*0.08, color="#e6e6e3", lw=0, label="returned-light weight (a.u.)")
Aw = {}
for lab, L, col, ls in curves:
    Ap = A_prime_vs_theta(L); Aw[lab] = np.sum(w_th*Ap)
    axc.plot(thetas, Ap, color=col, lw=2.0 if ls == "-" else 1.6, ls=ls, label=lab)
summary = "Angle-weighted round-trip loss ⟨A′⟩\n" + "\n".join(f"{k}: {v*100:.1f}%" for k, v in Aw.items())
axc.text(0.012, 0.97, summary, transform=axc.transAxes, ha="left", va="top", fontsize=8,
         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#d0d0cc", lw=0.6))
axc.set_xlim(0, 90); axc.set_ylim(0, 0.55)
axc.set_xlabel("Angle of incidence in substrate, θ (deg)"); axc.set_ylabel("Round-trip loss  A′(θ) = 1 − R_LED")
axc.set_title("(c)  Round-trip loss of the OLED stack seen by recycled light (EL-weighted, s/p average)", loc="left", fontsize=10)
axc.legend(frameon=False, fontsize=8, loc="upper right", ncol=1)
for s in ("top", "right"): axc.spines[s].set_visible(False)
axc.grid(axis="y", color="#e6e6e3", lw=0.6)

draw_profile(axd[0], stack("Al", K_HI), BLUE,   "(d) Al reflector, TCO k = 0.02")
draw_profile(axd[1], stack("Ag", K_HI), ORANGE, "(d) Ag reflector, TCO k = 0.02")
draw_profile(axd[2], stack("Ag", K_LO), ORANGE, "(d) Ag reflector, TCO k = 0.005")
axd[0].set_ylabel("Net flux S_z / returned flux")
fig.text(0.07, 0.015, "Draft, textbook n,k: substrate n = 1.77, TCO 50 nm, lossless organics 300 nm, metal 100 nm; "
         "Lambertian return proxy; EL Gaussian 530 nm / 60 nm FWHM.", fontsize=7.5, color="#52514e")
fig.savefig("fig2cd_mock.png", dpi=170); print("saved fig2cd_mock.png"); print({k: round(v*100, 2) for k, v in Aw.items()})

# ============================================================ SUPPLEMENTARY: DBR behind a TCO
fig2 = plt.figure(figsize=(11.0, 6.4), facecolor="white")
gs2 = fig2.add_gridspec(2, 2, hspace=0.55, wspace=0.3, left=0.07, right=0.98, top=0.93, bottom=0.10)
axa = fig2.add_subplot(gs2[0, :]); axb = fig2.add_subplot(gs2[1, 0]); axg = fig2.add_subplot(gs2[1, 1])
sc = [("DBR behind TCO, air-backed", stack("DBR", K_HI, "air"), AQUA, "-"),
      ("DBR behind TCO, glass-backed", stack("DBR", K_HI, "glass"), AQUA, "--"),
      ("Ag reflector (reference)", stack("Ag", K_HI), ORANGE, "-")]
axa.fill_between(thetas, 0, w_th/w_th.max()*0.1, color="#e6e6e3", lw=0, label="returned-light weight (a.u.)")
Aw2 = {}
for lab, L, col, ls in sc:
    Ap = A_prime_vs_theta(L); Aw2[lab] = np.sum(w_th*Ap); axa.plot(thetas, Ap, color=col, lw=2.0 if ls == "-" else 1.6, ls=ls, label=lab)
for x, lab in [(np.degrees(np.arcsin(1/N_SUB)), "TIR, air-backed"), (np.degrees(np.arcsin(1.5/N_SUB)), "TIR, glass-backed")]:
    axa.axvline(x, color="#9a9a96", lw=0.8, ls=":"); axa.text(x+0.6, 0.30, lab, rotation=90, va="top", ha="left", fontsize=7.5, color="#52514e")
axa.text(46, 0.77, "leaks into glass\n(no TIR below 58°)", fontsize=7.5, color="#52514e", ha="center", va="top")
axa.text(0.012, 0.97, "Angle-weighted ⟨A′⟩\n" + "\n".join(f"{k}: {v*100:.1f}%" for k, v in Aw2.items()),
         transform=axa.transAxes, ha="left", va="top", fontsize=8, bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#d0d0cc", lw=0.6))
axa.set_xlim(0, 90); axa.set_ylim(0, 0.8); axa.legend(frameon=False, fontsize=8, loc="upper right")
axa.set_xlabel("Angle of incidence in substrate, θ (deg)"); axa.set_ylabel("Round-trip loss  A′(θ)")
axa.set_title("(a)  Dielectric mirror behind a transparent electrode: the rear medium is part of the reflector", loc="left", fontsize=10)
for s in ("top", "right"): axa.spines[s].set_visible(False)
axa.grid(axis="y", color="#e6e6e3", lw=0.6)
draw_profile(axb, stack("DBR", K_HI, "air"), AQUA, "(b) DBR behind TCO, air-backed", ymax=0.30)
draw_profile(axg, stack("DBR", K_HI, "glass"), AQUA, "(c) DBR behind TCO, glass-backed", ymax=0.45)
axb.set_ylabel("Net flux S_z / returned flux")
fig2.text(0.07, 0.015, "Draft, textbook n,k: ZnS/LiF quarter-wave pairs at 540 nm; other parameters as in Fig. 2.", fontsize=7.5, color="#52514e")
fig2.savefig("figS_dbr_mock.png", dpi=170); print("saved figS_dbr_mock.png"); print({k: round(v*100, 2) for k, v in Aw2.items()})
