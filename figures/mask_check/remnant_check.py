"""Organic opening measured on the unified-scale drawings (22.0 px/mm from the slot: 12.7 x 8.7 mm, centred on the slot)
-> x 6.10-14.80, y 5.15-17.85 in the 21x23 pattern frame.  Tolerance (anode shift) with/without the strip trim and with/
without the two IZO remnants circled on the left (approx. x 5.9-6.6 / y 5.0-7.8 and x 5.7-6.4 / y 8.5-10.9).
A remnant counts as cathode-connected when cathode metal lands on it outside the organic area; it then leaks if any part
of it is under the organic."""
import numpy as np
from matplotlib.path import Path
W, H, CX = 21.0, 23.0, 10.45
R = lambda x0, y0, x1, y1: [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
mir = lambda p: [(2 * CX - x, y) for x, y in p]
both = lambda ps: ps + [mir(p) for p in ps]
ORG = (6.10, 5.15, 14.80, 17.85)
STRIP = [(7.95, 0), (12.95, 0), (12.95, 2.5), (13.94, 3.46), (13.94, 19.54), (12.95, 20.5), (12.95, 23), (7.95, 23),
         (7.95, 20.5), (6.96, 19.54), (6.96, 3.46), (7.95, 2.5)]
SLOT = R(8.97, 4.97, 11.93, 18.03)
TRIM = [R(8.5, 0, 12.4, 23), R(7.5, 5.0, 13.4, 18.0)]
UNTRIM_ETCHED = [STRIP]                                    # earlier plan: windows x 3.9-6.7, strip untouched
ISL_A = both([R(3.96, 0, 5.5, 3.3), R(1.2, 8.99, 3.9, 10.95), R(1.2, 12.05, 3.9, 14.01), R(3.96, 19.7, 5.5, 23)])
ISL_B = both([R(3.96, 0, 6.46, 3.3), R(1.2, 8.99, 3.9, 10.95), R(1.2, 12.05, 3.9, 14.01), R(3.96, 19.7, 6.46, 23)])
REM = [R(5.9, 5.0, 6.6, 7.8), R(5.7, 8.5, 6.4, 10.9)]
TL = [(3.96, 0), (6.91, 0), (6.91, 2.5), (5.97, 2.5), (5.97, 6.0), (9.98, 6.0), (9.98, 8.0), (3.96, 8.0)]
CAT = both([TL, R(1.2, 8.99, 9.98, 10.95), R(1.2, 12.05, 9.98, 14.01), [(x, H - y) for x, y in TL]])
PADS = [R(8.0, 0, 12.9, 2.5), R(8.0, 20.5, 12.9, 23)]
D = 0.025; P = 100
xs = np.arange(-P, int(W / D) + P) * D + D / 2; ys = np.arange(-P, int(H / D) + P) * D + D / 2
XX, YY = np.meshgrid(xs, ys); PTS = np.c_[XX.ravel(), YY.ravel()]
def ras(ps, holes=()):
    m = np.zeros(len(PTS), bool)
    for p in ps: m |= Path(p).contains_points(PTS)
    for p in holes: m &= ~Path(p).contains_points(PTS)
    return m.reshape(XX.shape)
CATM, PADM, ORGM = ras(CAT), ras(PADS), ras([R(*ORG)]); COUT = CATM & ~ORGM
def sh(m, k, l):
    o = np.zeros_like(m); h, w = m.shape
    o[max(0, l):h + min(0, l), max(0, k):w + min(0, k)] = m[max(0, -l):h - max(0, l), max(0, -k):w - max(0, k)]; return o
s = np.arange(-2, 2.0001, 0.05)
def tol(central, isl, rem):
    A, I = ras(central, [SLOT]), ras(isl); REMS = [ras([r]) for r in rem]
    T = np.zeros((len(s), len(s)), int)
    for i, dy in enumerate(s):
        for j, dx in enumerate(s):
            k, l = int(round(dx / D)), int(round(dy / D)); a, im = sh(A, k, l), sh(I, k, l)
            short = (COUT & a).any() or (PADM & im).any()
            leak = (im & ORGM).any()
            for r in REMS:
                rr = sh(r, k, l)
                if (COUT & rr).any() and (rr & ORGM).any(): leak = True
            T[i, j] = 2 if short else (1 if leak else 0)
    return T
def rad(T, lvl):
    rr = np.hypot(*np.meshgrid(s, s)); b = rr[T >= lvl]; return b.min() if b.size else np.inf
i0 = len(s) // 2
for lab, c, isl, rem in (('trimmed strip, remnants removed', TRIM, ISL_A, []), ('trimmed strip, remnants left', TRIM, ISL_A, REM),
                         ('untrimmed (windows 3.9-6.7), remnants removed', UNTRIM_ETCHED, ISL_B, []),
                         ('untrimmed, remnants left', UNTRIM_ETCHED, ISL_B, REM)):
    T = tol(c, isl, rem)
    print('%-48s at (0,0): %s | short-free radius %.2f mm | short+leak-free radius %.2f mm' % (lab, ['ok', 'LEAK', 'SHORT'][T[i0, i0]], rad(T, 2), rad(T, 1)))
