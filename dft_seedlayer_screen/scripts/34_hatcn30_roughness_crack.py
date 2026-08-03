"""Does a 30 nm HATCN underlayer (rough + cracked) hurt the transmittance of the
HATCN/Ag/cap electrode, and by how much?

Script 29 answered only the COHERENT question (cracks change the interference
condition) for a bare HATCN film and found ~0.1%p -- negligible. That answer does
not carry over once Ag sits on top, because a crack removes the ANCHORING SITES
under the Ag. Three loss channels are separated here:

  (1) coherent  : crack = locally missing HATCN -> different TMM stack.
                  Area-weighted incoherent average over cracked / intact regions.
  (2) scattering: rough interfaces scatter light out of the specular beam.
                  Scalar (Rayleigh-Rice) factor exp[-(2*pi*sigma*dn/lambda)^2]
                  per interface; for the metal, exp[-(4*pi*sigma*n/lambda)^2].
                  NOTE: an integrating-sphere measurement recovers most of this;
                  a display sees it as haze. Specular T is the pessimistic bound.
  (3) morphology: over a crack the Ag has no nitrile anchor -> dewets into islands.
                  Maxwell-Garnett effective medium (f = 0.5, 2x physical thickness,
                  same treatment as scripts 19-20) -> blue LSPR absorption.
                  THIS is the dominant channel.

Baseline: glass / HATCN(30) / Ag(9) / HATCN(40 cap).
"""
import numpy as np, os, json

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
wl = np.linspace(400, 700, 301)

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  520.9, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
n_ag = np.interp(wl, JC_WL, AG_N) + 1j * np.interp(wl, JC_WL, AG_K)

def bruggeman(eps_i, eps_h, f):
    """Symmetric EMA. Maxwell-Garnett is a dilute-inclusion theory and at f = 0.5
    puts a spurious pole in the denominator (it made the 'island' film opaque at
    650 nm in the first run). Bruggeman is the right choice at this fill factor.
    Root selected for a passive medium (Im n >= 0)."""
    b = (3 * f - 1) * eps_i + (2 - 3 * f) * eps_h
    disc = np.sqrt(b ** 2 + 8 * eps_i * eps_h)
    e1, e2 = (b + disc) / 4, (b - disc) / 4
    return np.where(np.sqrt(e1).imag >= 0, e1, e2)

NK = {
    "Ag": n_ag,
    "AgIsl": np.sqrt(bruggeman(n_ag ** 2, np.full_like(wl, 1.95 + 0j) ** 2, 0.5)),
    "HATCN": np.full_like(wl, 1.95 + 0j, dtype=complex),
    "air": np.ones_like(wl, dtype=complex),
    "glass": np.full_like(wl, 1.52 + 0j, dtype=complex),
}

def abeles(ns, ds, lam):
    ns = np.conj(ns)
    n0, nsub = ns[0], ns[-1]
    M = np.eye(2, dtype=complex)
    for j in range(1, len(ns) - 1):
        d = 2 * np.pi * ns[j] * ds[j] / lam
        M = M @ np.array([[np.cos(d), 1j * np.sin(d) / ns[j]],
                          [1j * ns[j] * np.sin(d), np.cos(d)]], dtype=complex)
    B, C = M @ np.array([1.0, nsub], dtype=complex)
    R = abs((n0 * B - C) / (n0 * B + C)) ** 2
    T = 4 * n0.real * nsub.real / abs(n0 * B + C) ** 2
    return T, R

def spec(stack, inc="glass", sup="air"):
    mats = [inc] + [m for m, _ in stack] + [sup]
    ds = [0.0] + [d for _, d in stack] + [0.0]
    T = np.zeros_like(wl); R = np.zeros_like(wl)
    for i, lam in enumerate(wl):
        ns = np.array([NK[m][i] for m in mats])
        T[i], R[i] = abeles(ns, ds, lam)
    return T, R, 1 - T - R

def Vlam(l): return 1.019 * np.exp(-285.4 * ((l / 1000.0) - 0.559) ** 2)
W = Vlam(wl); W = W / W.sum()
def atw(a, l): return float(np.interp(l, wl, a))
def phot(a): return float((a * W).sum())
Tg = spec([])[0]
def rel(T, l=None): return 100 * (atw(T, l) / atw(Tg, l) if l else phot(T) / phot(Tg))

WLS = [450, 510, 550, 650]
def line(lab, T, A):
    return (f"{lab:<38}" + "".join(f"{rel(T,l):>8.1f}" for l in WLS) + f"{rel(T):>9.1f}"
            + "   |" + "".join(f"{100*atw(A,l):>7.2f}" for l in WLS) + f"{100*phot(A):>8.2f}")

hdr = (f"{'scenario':<38}" + "".join(f"{l:>8d}" for l in WLS) + f"{'phot':>9}"
       + "   |" + "".join(f"{l:>7d}" for l in WLS) + f"{'phot':>8}")

BASE_STACK = [("HATCN", 30.0), ("Ag", 9.0), ("HATCN", 40.0)]
Tb, Rb, Ab = spec(BASE_STACK)

print("=" * 118)
print("RELATIVE TRANSMITTANCE (%)                            |  ABSORPTION (%)")
print(hdr)
print("-" * 118)
print(line("BASELINE  HATCN30/Ag9/cap (smooth)", Tb, Ab))

# ---------- channel 1: coherent effect of cracks (Ag stays continuous) ----------
print("\n-- (1) coherent only: crack = missing HATCN, Ag still continuous --")
Tc, Rc, Ac = spec([("Ag", 9.0), ("HATCN", 40.0)])          # crack region
for f in [0.05, 0.10, 0.20, 0.40]:
    T = (1 - f) * Tb + f * Tc
    A = (1 - f) * Ab + f * Ac
    print(line(f"  crack area {100*f:.0f}%", T, A))

# ---------- channel 2: roughness ----------
# 2a. specular (coherent) loss at a rough interface. The correct scalar factor for
#     the TRANSMITTED beam is exp[-(2*pi*sigma*|n1-n2|/lambda)^2]; the first run
#     wrongly used the REFLECTION TIS form exp[-(4*pi*sigma*n/lambda)^2], which
#     overstated the loss by ~3x in the exponent.
print("\n-- (2a) roughness -> loss of SPECULAR beam (scattered light still exits: haze,")
print("        an integrating sphere recovers it -- this is NOT absorption) --")
def scat(sig, n1, n2):
    return np.exp(-((2 * np.pi * sig * np.abs(n1 - n2)) / wl) ** 2)
for sig in [2.0, 5.0, 10.0]:
    s = scat(sig, NK["HATCN"], n_ag) * scat(sig, n_ag, NK["HATCN"])   # both Ag faces
    print(line(f"  sigma = {sig:.0f} nm RMS (specular only)", Tb * s, Ab))

# 2b. thickness non-uniformity: conformal Ag over a rough HATCN also varies in
#     local thickness. Incoherent average over a gaussian thickness distribution.
print("\n-- (2b) Ag thickness scatter from a rough underlayer (incoherent average) --")
for sig_t in [1.0, 2.0, 3.0]:
    ts = np.linspace(max(9 - 3 * sig_t, 0.5), 9 + 3 * sig_t, 25)
    w = np.exp(-0.5 * ((ts - 9.0) / sig_t) ** 2); w /= w.sum()
    T = np.zeros_like(wl); A = np.zeros_like(wl)
    for t, wt in zip(ts, w):
        Tt, Rt, At = spec([("HATCN", 30.0), ("Ag", t), ("HATCN", 40.0)])
        T += wt * Tt; A += wt * At
    print(line(f"  Ag 9 nm +/- {sig_t:.0f} nm (1 sigma)", T, A))

# ---------- channel 3: Ag dewets over the cracks ----------
print("\n-- (3) crack -> Ag loses anchor sites -> dewets to islands (MG, f=0.5) --")
Td, Rd, Ad = spec([("AgIsl", 18.0), ("HATCN", 40.0)])       # island Ag over crack
for f in [0.05, 0.10, 0.20, 0.40]:
    T = (1 - f) * Tb + f * Td
    A = (1 - f) * Ab + f * Ad
    print(line(f"  crack area {100*f:.0f}%  (Ag islands there)", T, A))

# ---------- combined worst case ----------
print("\n-- combined: 20% cracked + Ag islands there + sigma = 5 nm (specular) --")
f = 0.20
s = scat(5.0, NK["HATCN"], n_ag) * scat(5.0, n_ag, NK["HATCN"])
T = ((1 - f) * Tb + f * Td) * s
A = (1 - f) * Ab + f * Ad
print(line("  worst case (specular)", T, A))
print(line("  worst case (total, sphere)", (1 - f) * Tb + f * Td, A))

print("\n-- reference: fully dewetted (island) Ag everywhere --")
print(line("  100% island Ag", Td, Ad))
print("=" * 118)

json.dump({"note": "see stdout tables"}, open(os.path.join(RUNS, "hatcn30_crack.json"), "w"))
