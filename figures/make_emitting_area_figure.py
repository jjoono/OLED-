# -*- coding: utf-8 -*-
"""Emitting area of a lossy OLED versus one built on the proposed design rule.

Two isometric device stacks.  Nothing about the glow is chosen by eye: both the
size and the brightness of the emitting area come from the same round-trip loss
model as Fig.1(a) (make_roundtrip_figure.py), which supplies the palette too.

    rho  = (1 - eta) * R_metal * T_TCO^2        power kept per round trip
    Lam  = 1 / ln(1/rho)                        lateral spreading length
    tot  = eta / (1 - rho)                      total power that escapes

    glow radius    proportional to Lam   -> 2.32x
    glow amplitude proportional to tot   -> 1.84x
"""
import math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_roundtrip_figure as F            # palette, single source of truth

W, H = 1200, 492
FONT = "Helvetica Neue, Helvetica, Arial, Liberation Sans, sans-serif"
EX, EY = 330.0, 38.0                          # isometric basis, "east"
DX, DY = 118.0, -132.0                        # isometric basis, "depth"
OC, GLASS, TCO, ORG, MET = 14.0, 28.0, 10.0, 12.0, 18.0
BANDS = [(0.0, OC, "#dfeaf2", F.GLASS_L, "Outcoupling structure"),
         (OC, OC + GLASS, F.GLASS_F, F.GLASS_L, "Glass substrate"),
         (OC + GLASS, OC + GLASS + TCO, F.TCO_F, F.TCO_L, "TCO anode"),
         (OC + GLASS + TCO, OC + GLASS + TCO + ORG, F.ORG_F, F.ORG_L, "Organic layers"),
         (OC + GLASS + TCO + ORG, OC + GLASS + TCO + ORG + MET, F.MET_F, F.MET_L,
          "Metal cathode")]
THICK = BANDS[-1][1]
ORIGIN = ((180.0, 280.0), (700.0, 280.0))

def f(v):
    s = "%.3f" % v
    return s.rstrip("0").rstrip(".") if "." in s else s

def dark(h, k=0.84):
    h = h.lstrip("#")
    return "#%02x%02x%02x" % tuple(max(0, min(255, int(int(h[i:i+2], 16) * k)))
                                   for i in (0, 2, 4))

def P(O, u, v, t=0.0):
    return (O[0] + u * EX + v * DX, O[1] + u * EY + v * DY + t)

def poly(pts, fill, stroke=None, lw=1.0, op=1.0):
    s = '<polygon points="%s" fill="%s"' % (
        " ".join("%s,%s" % (f(x), f(y)) for x, y in pts), fill)
    if stroke:
        s += ' stroke="%s" stroke-width="%s" stroke-linejoin="round"' % (stroke, f(lw))
    if op != 1.0:
        s += ' opacity="%s"' % f(op)
    return s + "/>"

def ring(O, r, stroke, lw, op):
    pts = [P(O, .5 + r * math.cos(2 * math.pi * i / 72),
             .5 + r * math.sin(2 * math.pi * i / 72)) for i in range(73)]
    return ('<path d="M %s Z" fill="none" stroke="%s" stroke-width="%s" '
            'opacity="%s" stroke-linejoin="round"/>'
            % (" L ".join("%s,%s" % (f(x), f(y)) for x, y in pts), stroke, f(lw), f(op)))

def arrow(p0, p1, w, color, op=1.0):
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    L = math.hypot(dx, dy) or 1.0
    ux, uy = dx / L, dy / L
    hl, hw = max(7.0, w * 2.7), max(5.0, w * 1.9)
    q = (p1[0] - ux * hl, p1[1] - uy * hl)
    px, py = -uy, ux
    return ('<g opacity="%s"><line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
            'stroke-width="%s" stroke-linecap="butt"/><polygon points="%s,%s %s,%s %s,%s" '
            'fill="%s"/></g>'
            % (f(op), f(p0[0]), f(p0[1]), f(q[0]), f(q[1]), color, f(w), f(p1[0]), f(p1[1]),
               f(q[0] + px * hw), f(q[1] + py * hw), f(q[0] - px * hw), f(q[1] - py * hw), color))

def txt(x, y, s, size=14, col=F.INK, anchor="start", weight="400"):
    return ('<text x="%s" y="%s" text-anchor="%s" fill="%s" font-size="%s" '
            'font-weight="%s" font-family="%s">%s</text>'
            % (f(x), f(y), anchor, col, f(size), weight, FONT, s))

# ------------------------------------------------------- loss model -> geometry
def model(eta, r_met, t_tco):
    rho = (1 - eta) * r_met * t_tco ** 2
    return 1.0 / math.log(1.0 / rho), eta / (1.0 - rho)

LAM_C, TOT_C = model(.30, .72, .90)
LAM_D, TOT_D = model(.30, .98, .995)
LAM_REF = 0.30                                   # design-rule spreading, in panel widths
CASES = [dict(O=ORIGIN[0], lam=LAM_REF * LAM_C / LAM_D, amp=TOT_C / TOT_D, gid="gC",
              title="Conventional OLED",
              sub="strong absorption — light dies out before it can spread"),
         dict(O=ORIGIN[1], lam=LAM_REF, amp=1.0, gid="gD",
              title="Designed by the proposed design rule",
              sub="weak absorption — light spreads across the panel and escapes")]

o = []
A = o.append
A('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
  % (W, H, W, H))

# emissive glow: a circle on the panel surface, so a radial gradient carried into
# the isometric plane by the same basis that builds the slab
A("<defs>")
for c in CASES:
    O, lam, amp = c["O"], c["lam"], c["amp"]
    stops = []
    for i in range(10):
        t = i / 9.0
        r = 0.62 * t
        val = amp * math.exp(-r / lam)
        col = "#ffffff" if t < .12 else ("#74d4c2" if t < .38 else F.RAY)
        stops.append('<stop offset="%s" stop-color="%s" stop-opacity="%s"/>'
                     % (f(t), col, f(min(0.95, val * 0.9))))
    A('<radialGradient id="%s" gradientUnits="userSpaceOnUse" cx="0.5" cy="0.5" r="0.62" '
      'gradientTransform="matrix(%s,%s,%s,%s,%s,%s)">%s</radialGradient>'
      % (c["gid"], f(EX), f(EY), f(DX), f(DY), f(O[0]), f(O[1]), "".join(stops)))
A("</defs>")
A('<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H))
A('<title>Emitting area of a lossy OLED versus one designed with the proposed '
  'design rule</title>')

for c in CASES:
    O, lam, amp = c["O"], c["lam"], c["amp"]
    A('<g>')
    # --- slab: top face, then the two visible side faces, band by band
    A(poly([P(O, 0, 0), P(O, 1, 0), P(O, 1, 1), P(O, 0, 1)], BANDS[0][2],
           BANDS[0][3], 1.2))
    for t1, t2, fill, line, _ in BANDS:
        A(poly([P(O, 0, 0, t1), P(O, 1, 0, t1), P(O, 1, 0, t2), P(O, 0, 0, t2)],
               fill, line, 1.0))
    for t1, t2, fill, line, _ in BANDS:
        A(poly([P(O, 1, 0, t1), P(O, 1, 1, t1), P(O, 1, 1, t2), P(O, 1, 0, t2)],
               dark(fill), dark(line, .9), 1.0))
    A(poly([P(O, 0, 0), P(O, 1, 0), P(O, 1, 0, THICK), P(O, 0, 0, THICK)], "none",
           F.MUTED, 0.0, 0.0))

    # --- emitting area on the top surface
    A(poly([P(O, 0, 0), P(O, 1, 0), P(O, 1, 1), P(O, 0, 1)], "url(#%s)" % c["gid"]))
    for r in (0.15, 0.30, 0.45):
        op = min(0.85, 2.0 * amp * math.exp(-r / lam))
        if op > 0.05:
            A(ring(O, r, "#ffffff", 1.4, op * 0.8))
    # lateral spreading, drawn in the plane of the panel
    reach = 0.40 * lam / LAM_REF
    for ang in (0, 90, 180, 270):
        a = math.radians(ang)
        A(arrow(P(O, .5 + .09 * math.cos(a), .5 + .09 * math.sin(a)),
                P(O, .5 + reach * math.cos(a), .5 + reach * math.sin(a)),
                1.8 * amp + 0.9, "#00725f", 0.8))
    # --- outcoupled light
    for r, angs in ((0.0, [0]), (0.19, [45, 135, 225, 315]), (0.37, [45, 135, 225, 315])):
        for ang in angs:
            x = amp * math.exp(-r / lam)
            if x < 0.07:
                continue
            a = math.radians(ang)
            p0 = P(O, .5 + r * math.cos(a), .5 + r * math.sin(a))
            A(arrow(p0, (p0[0], p0[1] - (30 + 62 * x ** 0.5)), 7.5 * x ** 0.6, F.RAY))

    cx = O[0] + (EX + DX) / 2.0
    A(txt(cx, 30, c["title"], 19, F.INK, "middle", "700"))
    A(txt(cx, 53, c["sub"], 14.5, F.MUTED, "middle"))
    A('</g>')

# ------------------------------------------------------------------ layer keys
O = ORIGIN[0]
for (t1, t2, _, _, name), y_lab in zip(BANDS, (282, 312, 339, 358, 380)):
    y_tip = O[1] + (t1 + t2) / 2.0
    A(txt(168, y_lab + 5, name, 13.5, F.MUTED, "end"))
    A('<path d="M 174,%s L 182,%s" fill="none" stroke="%s" stroke-width="0.9"/>'
      % (f(y_lab), f(y_tip), F.LEADER))

# ------------------------------------------------------------------ legend
LG = 438.0
A(arrow((352, LG + 7), (352, LG - 15), 6.0, F.RAY))
A(txt(368, LG + 3, "outcoupled light", 13.5, F.MUTED))
A(arrow((536, LG - 4), (594, LG - 4), 3.2, "#00725f", 0.85))
A(txt(606, LG + 3, "lateral spreading inside the substrate", 13.5, F.MUTED))
A(txt(600, LG + 28, "emitting-area size \u221d spreading length (%.1f\u00d7)"
      "\u2003\u00b7\u2003 brightness \u221d total extracted light (%.1f\u00d7)"
      % (LAM_D / LAM_C, TOT_D / TOT_C), 13, F.LEADER, "middle"))
A('</svg>')

dest = os.path.join(os.path.dirname(os.path.abspath(__file__)), "emitting_area.svg")
open(dest, "w").write("\n".join(o))
print("wrote %s\n  spreading %.2fx   total light %.2fx" % (dest, LAM_D / LAM_C, TOT_D / TOT_C))
