"""k and sheet resistance are not independent in a TCO, and that settles it.

The objection is fair: IZO really is amorphous as deposited, really does not
need a crystallisation anneal, and really does have k around 0.006 in the
visible. The lab's own files say so. So on absorption alone a 50 nm IZO beats an
8 nm silver, and script 119 records that.

But the k a TCO shows and the sheet resistance it delivers come from the same
electrons. Free-carrier absorption in the Drude limit is

    eps2 = N e^3 / (eps0 m*^2 mu omega^3)     and     rho = 1 / (N e mu)

so eliminating N,

    eps2 = e^2 / (eps0 m*^2 mu^2 rho omega^3)

Absorption is inversely proportional to resistivity. A TCO is transparent
BECAUSE it is resistive. Quoting k = 0.006 without quoting the 5e-4 Ohm cm it
came with is quoting half a measurement.

Mobility is the only free variable that improves both, which is why every TCO
paper chases it, and why amorphous oxides -- whose s-orbital conduction band is
insensitive to disorder -- are interesting at all. But mu enters squared and
sits near 40-50 cm2/Vs for good IZO, so there is not much left to win.

Only the free-carrier part of k scales with resistivity. The measured IZO
spectrum rises steeply below 500 nm -- 0.0012 at 550 nm but 0.011 at 450 and
0.053 at 400 -- and that is the absorption edge of a 3.5 eV gap, fixed by the
material and indifferent to doping. The two are separated before scaling: the
red side follows lambda^3 as free carriers must (0.0021 at 650 against 0.0012 at
550, where lambda^3 predicts 0.0020), so the free-carrier term is anchored there
and the remainder at each wavelength is the edge.

A first version of this script scaled a single k value from 550 nm and
concluded, against its own table, that a TCO thick enough to match the silver
absorbs more. It does not; using one wavelength understated the blue end by an
order of magnitude and the conclusion was wrong in the other direction.
"""
import csv, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
m = __import__("118_zns_cap_compare")

RHO_IZO_RT = 5.0e-4          # Ohm cm, room-temperature sputtered IZO, literature
N_TCO = 2.0
TARGETS = [(50.0, "device-sane"), (100.0, "thick"), (200.0, "very thick")]
AG_RS, AG_A_AIR, AG_A_MATCH = 9.1, 12.56, 5.66


def split_k(name, anchor=550.0):
    """Separate free-carrier k (scales as lambda^3, and with 1/rho) from the edge."""
    l, n, k = m.load(name)
    kfc0 = m.interp(anchor, l, k)
    edge = {}
    for w in range(400, 801, 5):
        kfc = kfc0 * (w / anchor)**3
        edge[w] = max(m.interp(float(w), l, k) - kfc, 0.0)
    nfun = lambda w: m.interp(w, l, n)
    return kfc0, edge, nfun


def k_at(w, kfc0, edge, scale, anchor=550.0):
    """Total k at wavelength w when resistivity is divided by `scale`."""
    e = edge.get(int(round(w / 5) * 5), 0.0)
    return kfc0 * scale * (w / anchor)**3 + e


def hemis(stack, nout):
    num = den = 0.0
    for th in m.ANGLES:
        wt = math.sin(math.radians(th)) * math.cos(math.radians(th))
        for w in range(400, 701, 20):
            num += stack.A(float(w), float(th), nout) * wt
            den += wt
    return num / den


def main():
    kfc0, edge, nfun = split_k("IZO")
    print(f"calibration: lab's measured IZO, free-carrier k(550) = {kfc0:.4f},")
    print(f"             paired with rho = {RHO_IZO_RT:.1e} Ohm cm (RT sputtered)\n")
    print("k separated into free-carrier and band-edge parts:")
    print(f"{'nm':>5} {'k measured':>11} {'free carrier':>13} {'band edge':>11}")
    l, _, kk = m.load("IZO")
    for w in (400, 450, 500, 550, 600, 650, 700):
        km = m.interp(float(w), l, kk)
        kf = kfc0 * (w / 550.0)**3
        print(f"{w:>5} {km:>11.4f} {kf:>13.4f} {max(km-kf,0):>11.4f}")
    print("\nOnly the middle column moves with doping. The edge does not.\n")

    print("resistivity pushed down at fixed mobility:")
    print(f"{'rho':>10} {'k_fc(550)':>10} {'Rs 50nm':>9} {'Rs 100nm':>9} {'Rs 200nm':>9}")
    for rho in (5e-4, 3e-4, 2e-4, 1e-4, 5e-5):
        sc = RHO_IZO_RT / rho
        print(f"{rho:>10.1e} {kfc0*sc:>10.4f} {rho/50e-7:>9.1f} {rho/100e-7:>9.1f} "
              f"{rho/200e-7:>9.1f}")
    print(f"\nBest ITO ever is ~1e-4 and needs crystallisation above 200 C;")
    print(f"room-temperature IZO is 3-5e-4.\n")

    print("one-pass absorption, ZnS capped, at each operating point\n")
    print(f"{'electrode':<26} {'rho':>9} {'Rs':>7} {'air':>9} {'matched':>9}")
    print(f"{'HATCN 5 / Ag 8':<26} {'7.3e-06':>9} {AG_RS:>6.1f} "
          f"{AG_A_AIR:>8.2f}% {AG_A_MATCH:>8.2f}%")
    rows = []
    for d, _ in TARGETS + [(550.0, "equal Rs")]:
        for rho in (RHO_IZO_RT, 1e-4):
            sc = RHO_IZO_RT / rho
            rs = rho / (d * 1e-7)
            st = m.Stack(f"IZO {d:.0f}",
                         [(lambda w, s_=sc: complex(nfun(w), k_at(w, kfc0, edge, s_)), d)])
            for _, no in m.N_OUT:
                st.optimise(no)
            a1, a2 = hemis(st, 1.0), hemis(st, 1.8)
            rows.append((d, rho, rs, a1, a2))
            lbl = f"IZO {d:.0f} nm" + ("" if rho == RHO_IZO_RT else " (best-ITO rho)")
            print(f"{lbl:<26} {rho:>9.1e} {rs:>6.1f} {a1*100:>8.2f}% {a2*100:>8.2f}%")

    # where the two cross, read off the rows just computed
    rt = [r for r in rows if r[1] == RHO_IZO_RT]
    below = [r for r in rt if r[3] * 100 < AG_A_AIR]
    above = [r for r in rt if r[3] * 100 >= AG_A_AIR]
    if below and above:
        lo, hi = max(below, key=lambda r: r[0]), min(above, key=lambda r: r[0])
        print(f"\nThe two cross between {lo[0]:.0f} nm (Rs {lo[2]:.0f}, "
              f"{lo[3]*100:.2f}%) and {hi[0]:.0f} nm")
        print(f"(Rs {hi[2]:.0f}, {hi[3]*100:.2f}%), so around Rs = 15 Ohm/sq. Above that")
        print("sheet resistance the oxide is the better electrode on both counts;")
        print("below it, only the silver gets there at all.")
    print("\nWHAT THIS MEANS FOR THE PAPER. 'Ultrathin silver beats a")
    print("room-temperature TCO optically' is false and is withdrawn. The")
    print("defensible framing is that the two occupy different points on a")
    print("conduction-transparency trade: the oxide is transparent BECAUSE it is")
    print("resistive, by the same Drude physics that ties the silver's absorption")
    print("to its own resistivity. Which end of that trade a device wants depends")
    print("on panel size and busbar pitch, and the paper should say so rather")
    print("than claim a general win.")


if __name__ == "__main__":
    main()
