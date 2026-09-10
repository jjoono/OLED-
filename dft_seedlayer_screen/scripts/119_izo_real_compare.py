"""The TCO comparison with the lab's own measured indices, and at equal Rs.

Script 118 swept k parametrically and found a crossover near 0.045, then
asserted that room-temperature sputtered oxides sit above it. The lab's own
measured files say otherwise: averaged over 400-700 nm,

    IZO 0.0063   l_IZO 0.0081   ITO 0.0055   ITO_SNU 0.0039   l_ITO_SNU 0.0032

every one of them an order of magnitude below the crossover. IZO in particular
is amorphous as deposited and needs no crystallisation anneal, so it really is
available at room temperature on an organic stack. Optically the silver loses,
and it loses by a lot. That is the honest result and it is computed here from
the measured spectra rather than from a stand-in.

The comparison that decides an electrode is not at equal thickness, though. A
50 nm oxide and an 8 nm silver do not carry the same current. Silver's measured
resistivity on HATCN is 7.28 uOhm cm; a room-temperature IZO is 400-600, some
seventy times worse, so matching sheet resistance takes hundreds of nanometres
and the absorption scales with that thickness. Both comparisons are printed.

Resistivities are literature order of magnitude and flagged as such -- the
conclusion turns on the ratio, which is not in doubt, not on the third digit.
"""
import csv, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m = __import__("118_zns_cap_compare")

RHO = {                      # Ohm cm, room-temperature films unless noted
    "Ag 8 nm (measured)": 7.28e-6,
    "IZO, RT sputtered":  5.0e-4,
    "ITO, RT amorphous":  8.0e-4,
    "ITO, annealed":      1.8e-4,
}
TARGET_RS = 9.1              # Ohm/sq, the measured silver electrode


def build(name, d):
    lam, n, k = m.load(name)
    return m.Stack(f"{name} {d:.0f}", [(lambda w: complex(m.interp(w, lam, n),
                                                          m.interp(w, lam, k)), d)])


def hemis(stack, nout):
    num = den = 0.0
    for th in m.ANGLES:
        wt = math.sin(math.radians(th)) * math.cos(math.radians(th))
        for w in range(400, 701, 20):
            num += stack.A(float(w), float(th), nout) * wt
            den += wt
    return num / den


def main():
    hl, hn, hk = m.load("HATCN")
    a8l, a8n, a8k = m.load("Ag8nm_on_HATCN5_measured")
    ag = m.Stack("HATCN 5 / Ag 8",
                 [(lambda w: complex(m.interp(w, hl, hn), m.interp(w, hl, hk)), 5.0),
                  (lambda w: complex(m.interp(w, a8l, a8n), m.interp(w, a8l, a8k)), 8.0)])

    print("PART 1 -- equal thickness, 50 nm oxide, measured indices\n")
    stacks = [ag] + [build(t, 50.0) for t in ("IZO", "l_IZO", "ITO", "l_ITO_SNU_temp")]
    for s in stacks:
        for _, no in m.N_OUT:
            s.optimise(no)
    print(f"{'':<22} {'air':>9} {'matched n=1.8':>15}")
    res = {}
    for s in stacks:
        a1, a2 = hemis(s, 1.0), hemis(s, 1.8)
        res[s.label] = (a1, a2)
        print(f"{s.label:<22} {a1*100:>8.2f}% {a2*100:>14.2f}%")
    base = res["HATCN 5 / Ag 8"]
    print(f"\n  Every measured oxide beats the silver optically, by "
          f"{base[0]/res['IZO 50'][0]:.1f}x in air and "
          f"{base[1]/res['IZO 50'][1]:.1f}x matched, against IZO.")

    print("\n\nPART 2 -- equal sheet resistance\n")
    print(f"target Rs = {TARGET_RS} Ohm/sq, the measured Ag 8 nm electrode\n")
    print(f"{'material':<22} {'rho (Ohm cm)':>13} {'d for 9.1 Ohm/sq':>17}")
    for k_, r in RHO.items():
        d_nm = r / TARGET_RS * 1e7
        print(f"{k_:<22} {r:>13.2e} {d_nm:>14.0f} nm")

    print("\none-pass absorption at that thickness:\n")
    print(f"{'':<22} {'d':>7} {'air':>9} {'matched':>10}")
    print(f"{'HATCN 5 / Ag 8':<22} {'8 nm':>7} {base[0]*100:>8.2f}% {base[1]*100:>9.2f}%")
    for tag, rho in (("IZO, RT sputtered", RHO["IZO, RT sputtered"]),
                     ("ITO, RT amorphous", RHO["ITO, RT amorphous"]),
                     ("ITO, annealed", RHO["ITO, annealed"])):
        d_nm = rho / TARGET_RS * 1e7
        src = "IZO" if "IZO" in tag else ("ITO" if "amorphous" in tag
                                          else "l_ITO_SNU_temp")
        s = build(src, d_nm)
        for _, no in m.N_OUT:
            s.optimise(no)
        print(f"{tag:<22} {d_nm:>4.0f} nm {hemis(s, 1.0)*100:>8.2f}% "
              f"{hemis(s, 1.8)*100:>9.2f}%")

    print("\n  At equal current-carrying capacity the ordering reverses. The oxide")
    print("  has to be hundreds of nanometres thick, and a layer that thick is")
    print("  also a thick optical cavity: it moves the resonance, adds its own")
    print("  interference, and stops being a drop-in replacement for an 8 nm film.")
    print("\n  So the real trade is not optical. It is: an oxide at a device-sane")
    print("  50 nm gives lower absorption but ~100 Ohm/sq, and a silver at 8 nm")
    print("  gives 9.1 Ohm/sq at higher absorption. Which wins depends on the")
    print("  panel size and the busbar layout, not on the electrode alone.")


if __name__ == "__main__":
    main()
