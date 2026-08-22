"""One-pass absorption from 0 to 80 degrees, for five top-electrode options.

Normal incidence is the least interesting angle in a top-emitting OLED. Light
generated in the organic hits the electrode over the full hemisphere, everything
past the critical angle is trapped and comes back for another pass, and an MLA
turns those repeat passes into the emission that actually escapes. An electrode
is therefore judged by the whole curve, not by its 0 degree value.

Definition, identical for every candidate: light arrives from the organic
(n = 1.8) at angle theta, crosses the electrode stack, and leaves into air.
A(theta) = 1 - T - R, averaged over s and p, which is the loss per pass. The
capping thickness is optimised once per electrode at normal incidence, exactly
as a real design would be, and then held fixed across the sweep.

Beyond the critical angle -- 33.7 degrees out of n = 1.8 into air -- T is zero
and A = 1 - R outright. That is the regime the recycling lives in, and it is
where the electrodes separate.

TCO INDEX ASSIGNMENT IS UNCONFIRMED. The library carries several ITO and IZO
entries whose deposition history was never recorded with them. The mapping used
here is inferred from the optical constants themselves -- an as-deposited
amorphous film is less dense and more defective, so lower n and higher k -- and
the alternates are printed alongside so the choice can be checked rather than
trusted.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
m = importlib.import_module("91_nk_library_series")

L = 550.0
N_ORG, N_CPL = 1.80, 2.10
ANGLES = np.arange(0.0, 80.01, 1.0)

# label -> (list of (nk, thickness_nm)), electrode layers only
def nk_at(name, lam=L):
    a = np.loadtxt(os.path.join(m.NK, f"{name}.csv"), delimiter=",", skiprows=1)
    return complex(np.interp(lam, a[:, 0], a[:, 1]), np.interp(lam, a[:, 0], a[:, 2]))


def measured_ag(seed, d, sid):
    T = m.read_csv(os.path.join(m.RAW, f"{sid}T.csv"))[L]
    R = m.read_csv(os.path.join(m.RAW, f"{sid}R.csv"))[L]
    nk, _ = m.invert(L, T, R, seed, float(d))
    return nk


def tmm(n, d, lam, theta0):
    """(R, T) for s and p at incidence theta0 in medium n[0]. d[0], d[-1] unused.

    Fresnel coefficients and the transmitted-flux factor have to come from the
    same convention or the two disagree near the critical angle, where cos goes
    to zero and any mismatch is amplified without limit. Both are the standard
    amplitude form here.
    """
    n = np.array([complex(x) for x in n])
    s0 = n[0] * np.sin(theta0)
    cos = np.sqrt(1.0 - (s0 / n)**2)
    out = []
    for pol in ("s", "p"):
        M = np.eye(2, dtype=complex)
        for j in range(len(n) - 1):
            ni, nj, ci, cj = n[j], n[j + 1], cos[j], cos[j + 1]
            if pol == "s":
                r = (ni * ci - nj * cj) / (ni * ci + nj * cj)
                t = 2 * ni * ci / (ni * ci + nj * cj)
            else:
                r = (nj * ci - ni * cj) / (nj * ci + ni * cj)
                t = 2 * ni * ci / (nj * ci + ni * cj)
            I = np.array([[1, r], [r, 1]], dtype=complex) / t
            if j + 1 < len(n) - 1:
                dl = 2 * np.pi / lam * nj * cj * d[j + 1]
                I = I @ np.array([[np.exp(-1j * dl), 0], [0, np.exp(1j * dl)]])
            M = M @ I
        r = M[1, 0] / M[0, 0]
        t = 1.0 / M[0, 0]
        R = abs(r)**2
        if pol == "s":
            T = abs(t)**2 * (n[-1] * cos[-1]).real / (n[0] * cos[0]).real
        else:
            T = abs(t)**2 * (np.conj(n[-1]) * cos[-1]).real / (np.conj(n[0]) * cos[0]).real
        out.append((R, max(T, 0.0)))
    return out


def absorb(layers, d_cpl, theta):
    n = [N_ORG] + [x[0] for x in layers] + [N_CPL, 1.0]
    d = [0.0] + [x[1] for x in layers] + [d_cpl, 0.0]
    (Rs_, Ts_), (Rp, Tp) = tmm(n, d, L, np.deg2rad(theta))
    return 0.5 * ((1 - Rs_ - Ts_) + (1 - Rp - Tp))


def best_cpl(layers):
    best = (9e9, 0.0)
    for dc in np.arange(20.0, 121.0, 1.0):   # a device-realistic capping range
        a = absorb(layers, dc, 0.0)
        if a < best[0]:
            best = (a, dc)
    return best[1]


def main():
    ag_h = measured_ag("HATCN", 8, "2-2")
    ag_m = measured_ag("MoOx", 8, "2-10")

    stacks = [
        ("HATCN 5 / Ag 8",   [(nk_at("HATCN"), 5.0), (ag_h, 8.0)]),
        ("MoOx 5 / Ag 8",    [(nk_at("MoO3"), 5.0),  (ag_m, 8.0)]),
        ("ITO 70 (열처리 X)", [(nk_at("ITO"), 70.0)]),
        ("ITO 70 (열처리 O)", [(nk_at("l_ITO_SNU_temp"), 70.0)]),
        ("IZO 70",           [(nk_at("l_IZO"), 70.0)]),
    ]

    print(f"one-pass absorption at {L:.0f} nm, incident from organic n = {N_ORG}, "
          f"CPL n = {N_CPL}")
    print(f"critical angle into air: {np.rad2deg(np.arcsin(1.0/N_ORG)):.1f} deg\n")
    print("electrode optical constants used")
    for lbl, layers in stacks:
        top = layers[-1][0]
        print(f"  {lbl:<20} n = {top.real:6.3f}, k = {top.imag:6.4f}")
    print("\n  alternates not used (assignment unconfirmed):")
    for nm in ("ITO_SNU", "IZO"):
        v = nk_at(nm)
        print(f"    {nm:<16} n = {v.real:6.3f}, k = {v.imag:6.4f}")

    cpl = {lbl: best_cpl(layers) for lbl, layers in stacks}
    print("\noptimised capping thickness")
    for lbl, _ in stacks:
        print(f"  {lbl:<20} {cpl[lbl]:.0f} nm")

    curves = {lbl: np.array([absorb(layers, cpl[lbl], t) for t in ANGLES])
              for lbl, layers in stacks}

    print(f"\n{'angle':>6} " + " ".join(f"{lbl:>20}" for lbl, _ in stacks))
    print("-" * (7 + 21 * len(stacks)))
    for t in (0, 10, 20, 30, 34, 40, 50, 60, 70, 80):
        i = int(np.argmin(np.abs(ANGLES - t)))
        row = " ".join(f"{curves[lbl][i]*100:>19.2f}%" for lbl, _ in stacks)
        print(f"{t:>5}° {row}")

    # what the device actually experiences: sin-weighted over the hemisphere
    w = np.sin(np.deg2rad(ANGLES)) * np.cos(np.deg2rad(ANGLES))
    print(f"\n{'':>6} " + " ".join(f"{lbl:>20}" for lbl, _ in stacks))
    for name, sel in (("0-34° (escape cone)", ANGLES <= 33.75),
                      ("34-80° (trapped)", ANGLES > 33.75),
                      ("0-80° (all modes)", ANGLES <= 80.0)):
        vals = []
        for lbl, _ in stacks:
            v = np.sum(curves[lbl][sel] * w[sel]) / np.sum(w[sel])
            vals.append(f"{v*100:>19.2f}%")
        print(f"{name:>19} " + " ".join(vals))

    print("\nratio to HATCN/Ag, averaged over the trapped modes:")
    sel = ANGLES > 33.75
    base = np.sum(curves["HATCN 5 / Ag 8"][sel] * w[sel]) / np.sum(w[sel])
    for lbl, _ in stacks:
        v = np.sum(curves[lbl][sel] * w[sel]) / np.sum(w[sel])
        print(f"  {lbl:<20} {v/base:5.2f}x")

    out = os.path.join(BASE_OUT := os.path.join(m.BASE, "data"), "angular_onepass.csv")
    with open(out, "w") as f:
        f.write("angle_deg," + ",".join(lbl.replace(",", " ") for lbl, _ in stacks) + "\n")
        for i, t in enumerate(ANGLES):
            f.write(f"{t:.0f}," + ",".join(f"{curves[lbl][i]*100:.4f}"
                                           for lbl, _ in stacks) + "\n")
    print(f"\nwrote {os.path.relpath(out, m.BASE)}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    COL = {"HATCN 5 / Ag 8": "#2EC4B6", "MoOx 5 / Ag 8": "#FF6B35",
           "ITO 70 (열처리 X)": "#8899AA", "ITO 70 (열처리 O)": "#1E2761",
           "IZO 70": "#7B5EA7"}
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6),
                           gridspec_kw={"width_ratios": [2.05, 1]})
    EN = {"HATCN 5 / Ag 8": "HATCN 5 / Ag 8", "MoOx 5 / Ag 8": "MoOx 5 / Ag 8",
          "ITO 70 (열처리 X)": "ITO 70, as-deposited", "ITO 70 (열처리 O)": "ITO 70, annealed",
          "IZO 70": "IZO 70"}
    tc = np.rad2deg(np.arcsin(1.0 / N_ORG))
    ax[0].axvspan(tc, 80, color="#000000", alpha=0.045, lw=0)
    ax[0].axvline(tc, color="#555555", lw=0.9, ls="--")
    ax[0].text(tc + 1.2, 47, f"critical angle {tc:.1f}°\ntrapped, recycled by the MLA",
               fontsize=8.5, color="#555555", va="top")
    for lbl, _ in stacks:
        ax[0].plot(ANGLES, curves[lbl] * 100, color=COL[lbl], lw=2.2, label=EN[lbl])
    ax[0].set_xlabel("angle in the organic (deg)")
    ax[0].set_ylabel("one-pass absorption (%)")
    ax[0].set_xlim(0, 80); ax[0].set_ylim(0, 50)
    ax[0].legend(frameon=False, fontsize=9.5, loc="upper left")
    ax[0].set_title("Loss per pass through the top electrode", fontsize=11, loc="left")
    ax[0].grid(axis="y", color="#DDE3EC", lw=0.8)
    ax[0].set_axisbelow(True)
    for sp_ in ("top", "right"):
        ax[0].spines[sp_].set_visible(False)

    names = [EN[l] for l, _ in stacks]
    esc = [np.sum(curves[l][ANGLES <= tc] * w[ANGLES <= tc]) /
           np.sum(w[ANGLES <= tc]) * 100 for l, _ in stacks]
    trp = [np.sum(curves[l][ANGLES > tc] * w[ANGLES > tc]) /
           np.sum(w[ANGLES > tc]) * 100 for l, _ in stacks]
    y = np.arange(len(names))
    ax[1].barh(y - 0.2, esc, 0.38, color="#B8C4D4", label="escape cone 0-34°")
    ax[1].barh(y + 0.2, trp, 0.38, color=[COL[l] for l, _ in stacks])
    for i, (a, b) in enumerate(zip(esc, trp)):
        ax[1].text(a + 0.6, i - 0.2, f"{a:.1f}", va="center", fontsize=8.5, color="#555")
        ax[1].text(b + 0.6, i + 0.2, f"{b:.1f}", va="center", fontsize=8.5,
                   fontweight="bold")
    ax[1].set_yticks(y); ax[1].set_yticklabels(names, fontsize=9)
    ax[1].invert_yaxis()
    ax[1].set_xlabel("solid-angle weighted absorption (%)")
    ax[1].set_xlim(0, 44)
    # the trapped bars are coloured per electrode, so the legend key for them has
    # to be neutral rather than borrowing whichever colour came first
    from matplotlib.patches import Patch
    ax[1].legend(handles=[Patch(facecolor="#B8C4D4", label="escape cone 0-34°"),
                          Patch(facecolor="#6B7280", label="trapped 34-80°")],
                 frameon=False, fontsize=8.5, loc="lower right")
    ax[1].set_title("Where the loss actually happens", fontsize=11, loc="left")
    ax[1].grid(axis="x", color="#DDE3EC", lw=0.8); ax[1].set_axisbelow(True)
    for sp_ in ("top", "right", "left"):
        ax[1].spines[sp_].set_visible(False)
    fig.tight_layout()
    png = os.path.join(m.BASE, "data", "angular_onepass.png")
    fig.savefig(png, dpi=190)
    print(f"wrote {os.path.relpath(png, m.BASE)}")


if __name__ == "__main__":
    main()
