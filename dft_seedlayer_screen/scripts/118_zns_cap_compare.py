"""One-pass absorption with a ZnS cap, against TCO, resolved in wavelength and angle.

The electrode is the one actually made: HATCN 5 nm seed, silver at the measured
optical constants for 7 and 8 nm on that seed, capped with ZnS. Light arrives
from the organic (n = 1.8), crosses the electrode, and leaves into air, and
A = 1 - T - R averaged over s and p is the loss per pass.

ZnS is used because it is the highest index a thermally evaporated transparent
cap reaches, n = 2.39 at 550 nm against the 2.10 assumed for organic caps
earlier. Its dispersion is Debenham's Sellmeier for cubic ZnS; k is taken as
zero, the bandgap being 3.54 eV, so every ZnS number here is the ideal-cap limit
and a real evaporated film with residual absorption sits slightly above it.

The TCO comparison is parametric rather than taken from a measured file: n = 2.0
with k swept from 0.01 to 0.1 at 50 nm, because the deposition history of the
ITO and IZO entries in the library was never recorded and their k spans that
whole range. Sweeping says which part of the range a TCO must be in to beat the
metal, instead of asserting one value.

Cap thickness is optimised once per electrode at normal incidence and held fixed
across the angular sweep, as a real design would.
"""
import csv, math, os

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
NK = os.path.join(BASE, "data", "nk")
N_ORG = 1.80
ANGLES = list(range(0, 81))
TCO_D, TCO_N = 50.0, 2.0
TCO_K = [0.01, 0.03, 0.05, 0.1]


def load(name):
    lam, n, k = [], [], []
    with open(os.path.join(NK, name + ".csv")) as f:
        for row in list(csv.reader(f))[1:]:
            lam.append(float(row[0])); n.append(float(row[1])); k.append(float(row[2]))
    return lam, n, k


def interp(x, xs, ys):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    for i in range(1, len(xs)):
        if xs[i] >= x:
            t = (x - xs[i-1]) / (xs[i] - xs[i-1])
            return ys[i-1] + t * (ys[i] - ys[i-1])
    return ys[-1]


def zns_n(lam_nm):
    """Debenham Sellmeier for cubic ZnS, lambda in micron."""
    l2 = (lam_nm / 1000.0)**2
    return math.sqrt(8.393 + 0.14383 / (l2 - 0.2421**2) + 4430.99 / (l2 - 36.71**2))


def tmm(n, d, lam, th0):
    """(R, T) averaged over s and p; n complex, d in nm, th0 taken in medium 0."""
    s0 = n[0] * math.sin(th0)
    cos = [(1 - (s0 / x)**2)**0.5 for x in n]
    RT = []
    for pol in ("s", "p"):
        M = [[1 + 0j, 0j], [0j, 1 + 0j]]
        for j in range(len(n) - 1):
            ni, nj, ci, cj = n[j], n[j+1], cos[j], cos[j+1]
            if pol == "s":
                r = (ni*ci - nj*cj) / (ni*ci + nj*cj); t = 2*ni*ci / (ni*ci + nj*cj)
            else:
                r = (nj*ci - ni*cj) / (nj*ci + ni*cj); t = 2*ni*ci / (nj*ci + ni*cj)
            I = [[1/t, r/t], [r/t, 1/t]]
            if j + 1 < len(n) - 1:
                dl = 2*math.pi/lam * nj * cj * d[j+1]
                e, ei = complex(0, -1)**0, 0
                e = complex(math.e)**(-1j*dl) if False else None
                ph = -1j * dl
                e = complex(math.cos(ph.imag), math.sin(ph.imag)) * math.exp(ph.real)
                ph2 = 1j * dl
                ei = complex(math.cos(ph2.imag), math.sin(ph2.imag)) * math.exp(ph2.real)
                I = [[I[0][0]*e, I[0][1]*ei], [I[1][0]*e, I[1][1]*ei]]
            M = [[M[0][0]*I[0][0] + M[0][1]*I[1][0], M[0][0]*I[0][1] + M[0][1]*I[1][1]],
                 [M[1][0]*I[0][0] + M[1][1]*I[1][0], M[1][0]*I[0][1] + M[1][1]*I[1][1]]]
        r = M[1][0] / M[0][0]
        t = 1.0 / M[0][0]
        R = abs(r)**2
        if pol == "s":
            T = abs(t)**2 * (n[-1]*cos[-1]).real / (n[0]*cos[0]).real
        else:
            T = abs(t)**2 * (n[-1].conjugate()*cos[-1]).real / (n[0].conjugate()*cos[0]).real
        RT.append((R, max(T, 0.0)))
    return 0.5*(RT[0][0]+RT[1][0]), 0.5*(RT[0][1]+RT[1][1])


class Stack:
    def __init__(self, label, layers):
        self.label, self.layers, self.dcap = label, layers, 60.0

    def A(self, lam, th_deg, dcap=None):
        dc = self.dcap if dcap is None else dcap
        n = [complex(N_ORG, 0)] + [L(lam) for L, _ in self.layers] + \
            [complex(zns_n(lam), 0), complex(1.0, 0)]
        d = [0.0] + [t for _, t in self.layers] + [dc, 0.0]
        R, T = tmm(n, d, lam, math.radians(th_deg))
        return 1.0 - R - T

    def optimise(self, lam=550.0):
        best = min(((self.A(lam, 0.0, dc), dc) for dc in range(20, 141)),
                   key=lambda p: p[0])
        self.dcap = float(best[1])
        return self.dcap


def main():
    hl, hn, hk = load("HATCN")
    a7l, a7n, a7k = load("Ag7nm_on_HATCN5_measured")
    a8l, a8n, a8k = load("Ag8nm_on_HATCN5_measured")
    hatcn = lambda w: complex(interp(w, hl, hn), interp(w, hl, hk))
    ag7 = lambda w: complex(interp(w, a7l, a7n), interp(w, a7l, a7k))
    ag8 = lambda w: complex(interp(w, a8l, a8n), interp(w, a8l, a8k))

    stacks = [Stack("HATCN 5 / Ag 7", [(hatcn, 5.0), (ag7, 7.0)]),
              Stack("HATCN 5 / Ag 8", [(hatcn, 5.0), (ag8, 8.0)])]
    for k in TCO_K:
        stacks.append(Stack(f"TCO 50 k={k:.2f}",
                            [(lambda w, kk=k: complex(TCO_N, kk), TCO_D)]))
    for s in stacks:
        s.optimise()

    tc = math.degrees(math.asin(1 / N_ORG))
    print("one-pass A: organic n=1.8 -> electrode -> ZnS cap -> air")
    print(f"ZnS n = {zns_n(450):.3f} / {zns_n(550):.3f} / {zns_n(650):.3f} "
          f"at 450/550/650 nm, k = 0\n")
    print("optimised ZnS cap")
    for s in stacks:
        print(f"  {s.label:<18} {s.dcap:>5.0f} nm")

    print(f"\nA(lambda), normal incidence")
    print(f"{'nm':>5} " + " ".join(f"{s.label:>15}" for s in stacks))
    for w in range(400, 801, 50):
        print(f"{w:>5} " + " ".join(f"{s.A(w, 0.0)*100:>14.2f}%" for s in stacks))

    print(f"\nA(theta) at 550 nm   (critical angle {tc:.1f} deg)")
    print(f"{'deg':>5} " + " ".join(f"{s.label:>15}" for s in stacks))
    for th in (0, 20, 30, 34, 40, 50, 60, 70, 80):
        print(f"{th:>5} " + " ".join(f"{s.A(550.0, th)*100:>14.2f}%" for s in stacks))

    print(f"\nsolid-angle weighted, flat 400-700 nm")
    print(f"{'':<18} {'0-34 escape':>12} {'34-80 trapped':>14} {'0-80 all':>10}")
    rows = []
    for s in stacks:
        out = {}
        for name, sel in (("esc", lambda t: t <= tc), ("trap", lambda t: t > tc),
                          ("all", lambda t: True)):
            num = den = 0.0
            for th in ANGLES:
                if not sel(th):
                    continue
                wt = math.sin(math.radians(th)) * math.cos(math.radians(th))
                for w in range(400, 701, 20):
                    num += s.A(w, th) * wt
                    den += wt
            out[name] = num / den if den else 0.0
        rows.append((s.label, out))
        print(f"{s.label:<18} {out['esc']*100:>11.2f}% {out['trap']*100:>13.2f}% "
              f"{out['all']*100:>9.2f}%")

    ag = [r for r in rows if r[0].startswith("HATCN 5 / Ag 8")][0][1]["all"]
    print(f"\nk at which a 50 nm TCO matches HATCN 5 / Ag 8 ({ag*100:.2f}% over all modes):")
    lo, hi = 0.001, 0.5
    for _ in range(40):
        mid = (lo + hi) / 2
        s = Stack("probe", [(lambda w, kk=mid: complex(TCO_N, kk), TCO_D)])
        s.optimise()
        num = den = 0.0
        for th in ANGLES:
            wt = math.sin(math.radians(th)) * math.cos(math.radians(th))
            for w in range(400, 701, 20):
                num += s.A(w, th) * wt
                den += wt
        (lo, hi) = (mid, hi) if num / den < ag else (lo, mid)
    print(f"  k = {(lo+hi)/2:.3f}")
    out = os.path.join(BASE, "data", "zns_cap_onepass.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["quantity", "x"] + [s.label for s in stacks])
        for wl in range(400, 801, 5):
            w.writerow(["A_vs_wavelength_pct_at_0deg", wl] +
                       [f"{s.A(float(wl), 0.0)*100:.4f}" for s in stacks])
        for th in ANGLES:
            w.writerow(["A_vs_angle_pct_at_550nm", th] +
                       [f"{s.A(550.0, float(th))*100:.4f}" for s in stacks])
    print(f"\nwrote {os.path.relpath(out, BASE)}")
    print("ZnS k is taken as zero, so these are ideal-cap numbers.")


if __name__ == "__main__":
    main()
