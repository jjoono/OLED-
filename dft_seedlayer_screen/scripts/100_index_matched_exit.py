"""The same five electrodes with the total internal reflection removed.

Everything in script 99 assumes the light finally meets air, so a critical angle
sits at 33.7 degrees and everything past it is trapped. Raise the outcoupling
substrate to n = 1.8 -- matched to the organic -- and that boundary disappears:
no angle is trapped at the electrode any more, and the MLA on the substrate side
takes over the job of getting light out of n = 1.8 into air.

That changes what the electrode is being asked to do. With air outside, a
high-angle ray meets the electrode many times and each pass is a fresh chance to
be absorbed. With a matched substrate it crosses once and is gone, so the metal's
plasmon resonance is charged only once instead of on every bounce.

Two things are computed here, both against the air case for reference:

  A(theta)   loss on a single crossing, exit into n = 1.8 rather than 1.0
  effective  what the ray actually loses before it escapes, 1-(1-A)^S, with S
             the mean number of crossings

S is not 1 for the matched substrate. Nothing is trapped at the ELECTRODE any
more, but the substrate still has to give the light up to air, and whatever the
MLA fails to extract on a hit comes back down through the electrode. With a
per-hit escape probability q the mean crossing count is 1/q, so a realistic
MLA with an antireflection coat -- q around 0.75 -- costs about 1.3 crossings,
not 1. Three values of q are carried through rather than one, because that
number is the least certain thing here and the conclusion should be visible
across its plausible range.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
a99 = importlib.import_module("99_angular_seed_comparison")
m = a99.m

ANGLES = a99.ANGLES
N_MATCH = 1.80
PASSES = 5.9                    # mean passes before escape, air case, from the MLA analysis


def curves_for(stacks, n_out):
    cpl = {lbl: a99.best_cpl(layers, n_out) for lbl, layers in stacks}
    cur = {lbl: np.array([a99.absorb(layers, cpl[lbl], t, n_out) for t in ANGLES])
           for lbl, layers in stacks}
    return cpl, cur


def weighted(c, sel, w):
    return float(np.sum(c[sel] * w[sel]) / np.sum(w[sel]))


def main():
    stacks = a99.build_stacks()
    labels = [l for l, _ in stacks]
    w = np.sin(np.deg2rad(ANGLES)) * np.cos(np.deg2rad(ANGLES))
    tc = np.rad2deg(np.arcsin(1.0 / a99.N_ORG))

    cpl_air, cur_air = curves_for(stacks, 1.0)
    cpl_mat, cur_mat = curves_for(stacks, N_MATCH)

    print(f"one-pass absorption at {a99.L:.0f} nm, incident from organic "
          f"n = {a99.N_ORG}, CPL n = {a99.N_CPL}")
    print(f"  air exit      n = 1.00, critical angle {tc:.1f} deg")
    print(f"  matched exit  n = {N_MATCH:.2f}, no critical angle\n")

    print("optimised capping thickness")
    for l in labels:
        print(f"  {l:<20} air {cpl_air[l]:>3.0f} nm    matched {cpl_mat[l]:>3.0f} nm")

    print(f"\nA(theta) with the matched substrate")
    print(f"{'angle':>6} " + " ".join(f"{l:>20}" for l in labels))
    print("-" * (7 + 21 * len(labels)))
    for t in (0, 20, 34, 40, 50, 60, 70, 80):
        i = int(np.argmin(np.abs(ANGLES - t)))
        print(f"{t:>5}° " + " ".join(f"{cur_mat[l][i]*100:>19.2f}%" for l in labels))

    print(f"\nsingle-pass loss, solid-angle weighted over 0-80 deg")
    print(f"{'':>20} {'air exit':>10} {'matched':>10} {'change':>10}")
    sel = ANGLES <= 80.0
    for l in labels:
        va, vm = weighted(cur_air[l], sel, w), weighted(cur_mat[l], sel, w)
        print(f"{l:>20} {va*100:>9.2f}% {vm*100:>9.2f}% {(vm-va)*100:>+9.2f}%")

    print(f"\nloss before escape -- the quantity EQE sees")
    print(f"  air:     trapped past {tc:.1f} deg, about {PASSES:.1f} crossings there")
    print(f"  matched: nothing trapped at the electrode, but the MLA still has to")
    print(f"           clear n = {N_MATCH} into air; q = per-hit escape, S = 1/q")
    Q = [(1.00, "q = 1.00, ideal"), (0.75, "q = 0.75, MLA + AR"), (0.50, "q = 0.50, plain MLA")]
    hdr = f"{'':>20} {'air, recycled':>14}" + "".join(f"{lab:>21}" for _, lab in Q)
    print(hdr)
    rows = []
    for l in labels:
        eff_air = np.where(ANGLES > tc, 1 - (1 - cur_air[l])**PASSES, cur_air[l])
        va = weighted(eff_air, sel, w)
        vs = []
        for q, _ in Q:
            S = 1.0 / q
            vs.append(weighted(1 - (1 - cur_mat[l])**S, sel, w))
        rows.append((l, va, vs[1], vs))
        print(f"{l:>20} {va*100:>13.2f}%" + "".join(f"{v*100:>20.2f}%" for v in vs))

    h = [r for r in rows if r[0].startswith("HATCN")][0]
    mo = [r for r in rows if r[0].startswith("MoOx")][0]
    it = [r for r in rows if "열처리 O" in r[0]][0]
    print(f"\nAt q = 0.75, HATCN/Ag goes from {h[1]*100:.1f}% of the light lost at the "
          f"electrode to {h[2]*100:.1f}%.")
    print(f"MoOx/Ag goes from {mo[1]*100:.1f}% to {mo[2]*100:.1f}%.")
    print(f"The seed-layer gap, {abs(h[1]-mo[1])*100:.1f} %p with air, narrows to "
          f"{abs(h[2]-mo[2])*100:.1f} %p once the recycling is gone --")
    print("removing the trap is worth more than any seed, and it also shrinks the")
    print("prize for choosing between seeds.")
    print(f"The metal-to-TCO gap closes the same way: {(h[1]-it[1])*100:.1f} %p becomes "
          f"{(h[2]-it[2])*100:.1f} %p.")

    out = os.path.join(m.BASE, "data", "angular_onepass_matched.csv")
    with open(out, "w") as f:
        f.write("angle_deg," + ",".join(f"{l.replace(',',' ')}_air" for l in labels)
                + "," + ",".join(f"{l.replace(',',' ')}_n1.8" for l in labels) + "\n")
        for i, t in enumerate(ANGLES):
            f.write(f"{t:.0f}," + ",".join(f"{cur_air[l][i]*100:.4f}" for l in labels)
                    + "," + ",".join(f"{cur_mat[l][i]*100:.4f}" for l in labels) + "\n")
    print(f"\nwrote {os.path.relpath(out, m.BASE)}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
    COL = {"HATCN 5 / Ag 8": "#2EC4B6", "MoOx 5 / Ag 8": "#FF6B35",
           "ITO 70 (열처리 X)": "#8899AA", "ITO 70 (열처리 O)": "#1E2761",
           "IZO 70": "#7B5EA7"}
    EN = {"HATCN 5 / Ag 8": "HATCN 5 / Ag 8", "MoOx 5 / Ag 8": "MoOx 5 / Ag 8",
          "ITO 70 (열처리 X)": "ITO 70, as-deposited",
          "ITO 70 (열처리 O)": "ITO 70, annealed", "IZO 70": "IZO 70"}

    fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6),
                           gridspec_kw={"width_ratios": [1.5, 1.5, 1.15]})
    for k, (cur, ttl, shade) in enumerate(
            [(cur_air, "Air outside  (n = 1.0)", True),
             (cur_mat, f"Matched substrate  (n = {N_MATCH})", False)]):
        if shade:
            ax[k].axvspan(tc, 80, color="#000000", alpha=0.045, lw=0)
            ax[k].axvline(tc, color="#555555", lw=0.9, ls="--")
            ax[k].text(tc + 1.2, 47, f"trapped past {tc:.1f}°", fontsize=8.5,
                       color="#555555", va="top")
        else:
            ax[k].text(3, 47, "no trapped region", fontsize=8.5, color="#555555",
                       va="top")
        for l in labels:
            ax[k].plot(ANGLES, cur[l] * 100, color=COL[l], lw=2.2, label=EN[l])
        ax[k].set_xlabel("angle in the organic (deg)")
        ax[k].set_xlim(0, 80); ax[k].set_ylim(0, 50)
        ax[k].set_title(ttl, fontsize=11, loc="left")
        ax[k].grid(axis="y", color="#DDE3EC", lw=0.8); ax[k].set_axisbelow(True)
        for sp_ in ("top", "right"):
            ax[k].spines[sp_].set_visible(False)
    ax[0].set_ylabel("one-pass absorption (%)")
    ax[0].legend(frameon=False, fontsize=9.5, loc="upper left")

    y = np.arange(len(labels))
    va = [r[1] * 100 for r in rows]
    vm = [r[2] * 100 for r in rows]
    ax[2].barh(y - 0.2, va, 0.38, color="#B8C4D4")
    ax[2].barh(y + 0.2, vm, 0.38, color=[COL[l] for l in labels])
    for i, (p, q) in enumerate(zip(va, vm)):
        ax[2].text(p + 1.0, i - 0.2, f"{p:.1f}", va="center", fontsize=8.5, color="#555")
        ax[2].text(q + 1.0, i + 0.2, f"{q:.1f}", va="center", fontsize=8.5,
                   fontweight="bold")
    ax[2].set_yticks(y); ax[2].set_yticklabels([EN[l] for l in labels], fontsize=9)
    ax[2].invert_yaxis()
    ax[2].set_xlabel("light lost at the electrode before escape (%)")
    ax[2].set_xlim(0, max(va) * 1.22)
    ax[2].legend(handles=[Patch(facecolor="#B8C4D4", label=f"air, ~{PASSES:.1f} crossings"),
                          Patch(facecolor="#6B7280", label="matched, q = 0.75")],
                 frameon=False, fontsize=8.5, loc="lower right")
    ax[2].set_title("What the EQE actually pays", fontsize=11, loc="left")
    ax[2].grid(axis="x", color="#DDE3EC", lw=0.8); ax[2].set_axisbelow(True)
    for sp_ in ("top", "right", "left"):
        ax[2].spines[sp_].set_visible(False)
    fig.tight_layout()
    png = os.path.join(m.BASE, "data", "angular_matched_substrate.png")
    fig.savefig(png, dpi=185)
    print(f"wrote {os.path.relpath(png, m.BASE)}")


if __name__ == "__main__":
    main()
