"""What is the organic seed actually worth optically? TMM against metal seeds.

THE COMPARISON THAT MATTERS. Measurement gives HATCN(4 nm)/Ag(3 nm) = 100 Ohm/sq,
and the literature gives Cu-seeded Ag(3 nm) = 66 Ohm/sq. On sheet resistance the
metal seed wins. The case for an organic seed is therefore not electrical -- it is
that the metal seed buys its percolation with absorption, and this quantifies how
much.

WHY Cu IS HANDLED PARAMETRICALLY. The project has verified Johnson & Christy data
for Ag and Au and Palik/Rakic for Al, but not for Cu, and this session could not
retrieve a Cu table it could check (refractiveindex.info is blocked by the egress
proxy and the search results did not carry the numbers). Rather than type Cu's n
and k from memory into a calculation whose whole purpose is a quantitative
comparison, the generic-metal sweep below covers it: a 1 nm seed with k anywhere
in the 2-4 range that visible-wavelength Cu occupies is bracketed, and the answer
does not hinge on which value is right. Substitute a real Cu table when one is at
hand -- METALS below is where it goes.

CONFIGURATIONS. Both the bare stack that a spectrophotometer sees and the
in-device stack with an organic superstrate, because the seed sits at a different
field position in each and its absorption is not the same in the two.

HATCN's k is taken as 2e-4, the bound computed in scripts/47 from TD-DFT rather
than the k = 0 that scripts/32 assumed. That assumption is why the earlier
"+0.00 %p" figure was circular, and it is not repeated here.
"""
import os
import numpy as np

WL = np.arange(400.0, 701.0, 2.0)

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
AU_N = np.array([1.66, 1.70, 1.72, 1.76, 1.70, 1.13, 0.61, 0.37, 0.24, 0.17, 0.15, 0.14])
AU_K = np.array([1.96, 1.92, 1.90, 1.92, 2.00, 2.21, 2.50, 2.72, 3.09, 3.42, 3.86, 4.28])
AL_WL = np.array([400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0])
AL_N = np.array([0.49, 0.62, 0.77, 0.96, 1.20, 1.47, 1.83])
AL_K = np.array([4.86, 5.47, 6.08, 6.69, 7.26, 7.79, 8.31])

N_GLASS, N_AIR, N_ORG = 1.52, 1.0, 1.80
N_HATCN, K_HATCN = 1.95, 2.0e-4          # k from scripts/47, not assumed zero

AG_T = (3.0, 5.0)
SEED_METAL_T = 1.0                        # typical seed thickness in the literature
HATCN_T = 4.0                             # the measured stack
K_SWEEP = (1.0, 2.0, 3.0, 4.0)            # brackets Cu across the visible


def nk(name):
    if name == "Ag":
        return np.conj(np.interp(WL, JC_WL, AG_N) + 1j * np.interp(WL, JC_WL, AG_K))
    if name == "Au":
        return np.conj(np.interp(WL, JC_WL, AU_N) + 1j * np.interp(WL, JC_WL, AU_K))
    if name == "Al":
        return np.conj(np.interp(WL, AL_WL, AL_N) + 1j * np.interp(WL, AL_WL, AL_K))
    if name == "HATCN":
        return np.conj(np.full_like(WL, N_HATCN + 1j * K_HATCN, dtype=complex))
    raise KeyError(name)


def generic_metal(k, n=1.0):
    """A featureless metal seed: constant n and k across the visible."""
    return np.conj(np.full_like(WL, n + 1j * k, dtype=complex))


def tmm_layers(ns, ds, lam):
    """Abeles with layer-resolved absorption from the Poynting difference."""
    M = np.eye(2, dtype=complex)
    Ms = []
    for N, d in zip(ns[1:-1], ds[1:-1]):
        delta = 2 * np.pi / lam * N * d
        c, s = np.cos(delta), np.sin(delta)
        m = np.array([[c, 1j * s / N], [1j * N * s, c]])
        Ms.append(m)
        M = M @ m
    n0, nsub = ns[0], ns[-1]
    B, C = M @ np.array([1.0, nsub], dtype=complex)
    T = float(4 * n0.real * nsub.real / abs(n0 * B + C) ** 2)
    R = float(abs((n0 * B - C) / (n0 * B + C)) ** 2)

    # net Poynting flux at each interface, walking inward
    E, H = 1.0, n0 * (1 - (n0 * B - C) / (n0 * B + C)) / (1 + (n0 * B - C) / (n0 * B + C))
    Bv, Cv = B, C
    flux = []
    acc = np.eye(2, dtype=complex)
    for m in Ms:
        bc = acc @ np.array([1.0, nsub], dtype=complex)
        flux.append(float((np.conj(bc[0]) * bc[1]).real))
        acc = acc @ np.linalg.inv(m) if False else acc
    # simpler and numerically safe: peel matrices from the substrate side
    flux = []
    acc = np.array([1.0, nsub], dtype=complex)
    flux.append(float((np.conj(acc[0]) * acc[1]).real))
    for m in reversed(Ms):
        acc = m @ acc
        flux.append(float((np.conj(acc[0]) * acc[1]).real))
    flux = flux[::-1]                       # front -> back
    norm = abs(n0 * B + C) ** 2 / 4.0
    A = [(flux[i] - flux[i + 1]) / norm for i in range(len(Ms))]
    return T, R, A


def photopic(a):
    v = 1.019 * np.exp(-285.4 * ((WL / 1000.0) - 0.559) ** 2)
    return float((a * v).sum() / v.sum())


def run(stack, superstrate="air"):
    """stack = [(nk_array, d_nm), ...] between glass and the superstrate."""
    n_sup = N_AIR if superstrate == "air" else N_ORG
    T = np.zeros_like(WL)
    A = np.zeros((len(WL), len(stack)))
    for i, lam in enumerate(WL):
        ns = [complex(N_GLASS)] + [s[0][i] for s in stack] + [complex(n_sup)]
        ds = [0.0] + [s[1] for s in stack] + [0.0]
        t, r, a = tmm_layers(ns, ds, lam)
        T[i] = t
        A[i] = a
    return photopic(T), [photopic(A[:, j]) for j in range(len(stack))]


def main():
    ag = nk("Ag")
    for sup, label in (("air", "measurable (glass/seed/Ag/air)"),
                       ("org", "in device (glass/seed/Ag/organic)")):
        print("\n" + "=" * 74)
        print(label)
        print("=" * 74)
        for d_ag in AG_T:
            print(f"\n  Ag {d_ag:.0f} nm")
            print(f"    {'seed':<22}{'T (%)':>9}{'A_seed (%)':>12}"
                  f"{'A_Ag (%)':>10}{'dT vs HATCN':>14}")
            rows = {}
            # organic reference first
            T, A = run([(nk("HATCN"), HATCN_T), (ag, d_ag)], sup)
            rows["HATCN 4 nm"] = (T, A[0], A[1])
            for m in ("Al", "Au"):
                T2, A2 = run([(nk(m), SEED_METAL_T), (ag, d_ag)], sup)
                rows[f"{m} {SEED_METAL_T:.0f} nm"] = (T2, A2[0], A2[1])
            for k in K_SWEEP:
                T2, A2 = run([(generic_metal(k), SEED_METAL_T), (ag, d_ag)], sup)
                rows[f"metal k={k:.0f}, 1 nm"] = (T2, A2[0], A2[1])
            T0, A0 = run([(ag, d_ag)], sup)
            rows["none"] = (T0, 0.0, A0[0])

            tref = rows["HATCN 4 nm"][0]
            for name, (T, As, Aa) in rows.items():
                dt = 100 * (T - tref)
                mark = "  (reference)" if name.startswith("HATCN") else f"{dt:>+9.2f} %p"
                print(f"    {name:<22}{100*T:>9.2f}{100*As:>12.2f}{100*Aa:>10.2f}"
                      f"{mark:>14}")

    print("\n" + "=" * 74)
    print("READING THIS")
    print("=" * 74)
    print("""
  A_seed is the seed layer's OWN absorption, separated from the silver's. That
  column is the whole argument: it is what an organic seed does not pay and a
  metal seed does, at the same Ag thickness and therefore at comparable
  percolation.

  Cu is not tabulated here and the k sweep stands in for it. Across the visible
  Cu's k runs roughly 2-4, so the k=2 and k=4 rows bracket it; read the gain as
  the range between them rather than as a single number until a checked Cu table
  is substituted.

  HATCN's own absorption is not zero here -- it uses k = 2e-4 from scripts/47.
  It is small, but printing it rather than assuming it away is the point: the
  earlier "+0.00 %p" claim came from feeding k = 0 in and reading it back out.

  What this does NOT include: interfacial charge-transfer absorption at the
  HATCN/Ag junction, which is not a bulk k and which scripts/48 bounds
  separately at 0.02-0.9 % per interface depending on the CT oscillator
  strength. That is a real loss channel for the organic seed and it is not in
  these numbers.""")


if __name__ == "__main__":
    main()
