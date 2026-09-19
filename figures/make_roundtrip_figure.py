# -*- coding: utf-8 -*-
"""Fig.1(a): light recycling in an OLED with an external outcoupling structure.
Left  - conventional device: strong ohmic / TCO absorption, beam dies out.
Right - device from the proposed design rule: absorption suppressed, beam survives.
"""
import math, os

W, H = 2120, 442
FONT = "Helvetica Neue, Helvetica, Arial, Liberation Sans, sans-serif"

# ---------------------------------------------------------------- stack (y)
Y_CREST, Y_TOP = 145.0, 175.0      # lens apex / lens base = glass top surface
Y_TCO, Y_ORG, Y_MET, Y_BOT = 261.0, 277.0, 295.0, 335.0
Y_EMIT = 286.0                     # emitting plane, middle of the organic stack
LENS_R, PITCH = 30.0, 60.0
PW = 900.0                         # panel width = 15 lenses
PANEL_X = (200.0, 1160.0)

# ---------------------------------------------------------------- ray path
TAN = 1.0                          # 45 deg in the glass  (critical angle 41.8 deg)
HIT_TOP   = [150.0, 390.0, 630.0, 870.0]      # all of them lens centres
HIT_METAL = [270.0, 510.0, 750.0]
X_EMIT    = HIT_TOP[0] - (Y_EMIT - Y_TOP) * TAN

# ---- palette: Nature Publishing Group accents (ggsci "npg") over desaturated
# ---- structural tints.  Two saturated hues only: light vs. loss.
RAY       = "#00a087"   # NPG teal-green  - optical power
LOSS      = "#e64b35"   # NPG coral red   - absorption
BURST     = "#fbdcd6"   # pale tint of LOSS, fill of the absorption star
EMIT      = "#e8901f"   # gold            - exciton emission
EMIT_EDGE = "#6f4308"
GLASS_F, GLASS_L = "#eaf1f6", "#6e93ae"
TCO_F,   TCO_L   = "#cbe1ee", "#6e93ae"
ORG_F,   ORG_L   = "#fdf2e2", "#c99a55"
MET_F,   MET_L   = "#c3c7cb", "#868c92"
INK, MUTED, LEADER = "#1a1a1a", "#55585c", "#9aa0a6"
LW_SLAB, LW_LAYER, LW_STAR = 1.3, 1.05, 1.0
W_RAY = 9.5             # line width of the full-power ray
GREEN, BLUE, AMBER = RAY, LOSS, EMIT          # kept for backward compatibility

def f(v):
    s = "%.2f" % v
    return s.rstrip("0").rstrip(".") if "." in s else s

out = []
A = out.append

# ---------------------------------------------------------------- primitives
def ray(p0, p1, w, color=GREEN, op=1.0, head=True, trim0=0.0, trim1=0.0):
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    L = math.hypot(dx, dy) or 1.0
    ux, uy = dx / L, dy / L
    if trim0:
        p0 = (p0[0] + ux * trim0, p0[1] + uy * trim0)
    if trim1:
        p1 = (p1[0] - ux * trim1, p1[1] - uy * trim1)
    hl = max(7.0, w * 2.9)
    hw = max(5.0, w * 2.05)
    q = (p1[0] - ux * hl, p1[1] - uy * hl) if head else p1
    px, py = -uy, ux
    s = '<g opacity="%s">' % f(op)
    s += ('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" stroke-width="%s" '
          'stroke-linecap="butt"/>' % (f(p0[0]), f(p0[1]), f(q[0]), f(q[1]), color, f(w)))
    if head:
        s += ('<polygon points="%s,%s %s,%s %s,%s" fill="%s"/>'
              % (f(p1[0]), f(p1[1]), f(q[0] + px * hw), f(q[1] + py * hw),
                 f(q[0] - px * hw), f(q[1] - py * hw), color))
    return s + '</g>'

def star(cx, cy, r, fill, stroke, n=11, inner=0.44, op=1.0, sw=1.0):
    p = []
    for i in range(2 * n):
        rad = r if i % 2 == 0 else r * inner
        a = math.pi * i / n - math.pi / 2
        p.append((cx + rad * math.cos(a), cy + rad * math.sin(a)))
    return ('<polygon points="%s" fill="%s" stroke="%s" stroke-width="%s" '
            'stroke-linejoin="round" opacity="%s"/>'
            % (" ".join("%s,%s" % (f(x), f(y)) for x, y in p), fill, stroke, f(sw), f(op)))

def wave(x0, y, dx, amp, periods, color, sw, op, head=True):
    n, pt = 64, []
    for i in range(n + 1):
        t = i / float(n)
        pt.append((x0 + dx * t, y + amp * math.sin(2 * math.pi * periods * t)))
    d = "M " + " L ".join("%s,%s" % (f(x), f(y_)) for x, y_ in pt)
    s = ('<g opacity="%s"><path d="%s" fill="none" stroke="%s" stroke-width="%s" '
         'stroke-linecap="round" stroke-linejoin="round"/>' % (f(op), d, color, f(sw)))
    if head:
        sgn = 1.0 if dx > 0 else -1.0
        tip = (x0 + dx + sgn * 5.0, y)
        hw = max(3.4, sw * 2.0)
        s += ('<polygon points="%s,%s %s,%s %s,%s" fill="%s"/>'
              % (f(tip[0]), f(tip[1]), f(tip[0] - sgn * 8), f(y - hw),
                 f(tip[0] - sgn * 8), f(y + hw), color))
    return s + '</g>'

def txt(x, y, s, size=15, fill=INK, anchor="start", weight="400"):
    return ('<text x="%s" y="%s" text-anchor="%s" fill="%s" font-size="%s" '
            'font-weight="%s" font-family="%s">%s</text>'
            % (f(x), f(y), anchor, fill, f(size), weight, FONT, s))

# ---------------------------------------------------------------- one panel
def panel(x0, eta, r_met, t_tco, title, sub):
    g = ['<g id="panel">']
    # --- glass slab, then one closed semicircle per micro-lens
    g.append('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="%s" '
             'stroke-width="%s"/>' % (f(x0), f(Y_TOP), f(PW), f(Y_TCO - Y_TOP),
                                      GLASS_F, GLASS_L, f(LW_SLAB)))
    for k in range(int(PW / PITCH)):
        cx = x0 + PITCH * k + LENS_R
        g.append('<path d="M %s,%s A %s,%s 0 0 1 %s,%s Z" fill="%s" stroke="%s" '
                 'stroke-width="%s" stroke-linejoin="round"/>'
                 % (f(cx - LENS_R), f(Y_TOP), f(LENS_R), f(LENS_R), f(cx + LENS_R),
                    f(Y_TOP), GLASS_F, GLASS_L, f(LW_SLAB)))
    # --- device layers
    for y1, y2, fill, stroke in ((Y_TCO, Y_ORG, TCO_F, TCO_L),
                                 (Y_ORG, Y_MET, ORG_F, ORG_L),
                                 (Y_MET, Y_BOT, MET_F, MET_L)):
        g.append('<rect x="%s" y="%s" width="%s" height="%s" fill="%s" stroke="%s" '
                 'stroke-width="%s"/>' % (f(x0), f(y1), f(PW), f(y2 - y1), fill, stroke,
                                          f(LW_LAYER)))

    # --- radiometry along the zig-zag
    I, seg, esc, loss = 1.0, [], [], []
    for k in range(4):
        seg.append(I); esc.append(I * eta); I *= (1 - eta)
        if k < 3:
            I *= t_tco; loss.append(I * (1 - r_met)); I *= r_met; I *= t_tco
    wid = [W_RAY * v ** 0.75 for v in seg]

    # --- TCO absorption ticks, flanking every metal bounce
    tl = 1.0 - t_tco
    tw, ta = 10.0 + 180.0 * tl, 1.0 + 22.0 * tl
    for xm in HIT_METAL:
        for sgn in (-1, 1):
            g.append(wave(x0 + xm + sgn * 24 - tw / 2, Y_TCO + 8, tw, ta, 2, BLUE,
                          0.9 + 9.0 * tl, 0.25 + 6.5 * tl, head=False))

    # --- the ray itself
    kept = (1.0 - eta) ** 0.75           # power still travelling after a lens
    g.append(ray((x0 + X_EMIT, Y_EMIT), (x0 + HIT_TOP[0], Y_TOP), wid[0], trim0=10))
    for k in range(3):
        g.append(ray((x0 + HIT_TOP[k], Y_TOP), (x0 + HIT_METAL[k], Y_MET),
                     wid[k] * kept, trim0=13, trim1=9))
        g.append(ray((x0 + HIT_METAL[k], Y_MET), (x0 + HIT_TOP[k + 1], Y_TOP),
                     wid[k + 1], trim0=9))
    # tail: what is still bouncing after the last lens, fading out of frame
    g.append('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="url(#fade)" '
             'stroke-width="%s"/>' % (f(x0 + HIT_TOP[3] + 9), f(Y_TOP + 9), f(x0 + PW),
                                      f(Y_TOP + PW - HIT_TOP[3]), f(wid[3] * kept)))

    # --- escaping fans
    for k, xl in enumerate(HIT_TOP):
        fw = max(1.0, 0.55 * wid[k])
        for phi, ln in ((-30, 68), (-5, 79), (22, 70)):
            a = math.radians(phi)
            o = (x0 + xl + LENS_R * math.sin(a), Y_TOP - LENS_R * math.cos(a))
            b = math.radians(phi * 1.45)
            g.append(ray(o, (o[0] + ln * math.sin(b), o[1] - ln * math.cos(b)), fw))

    # --- absorption at the metal mirror
    for k, xm in enumerate(HIT_METAL):
        v = loss[k]
        r = 3.0 + 30.0 * v ** 0.6
        op = 0.5 + 0.5 * min(1.0, v / 0.18)
        ln = 10.0 + 200.0 * v
        sw = 1.2 + 6.0 * v ** 0.6
        for sgn in (-1, 1):
            g.append(wave(x0 + xm + sgn * (r - 1), Y_MET - 6, sgn * ln,
                          2.0 + 14 * v, 2.5, BLUE, sw, op))
        g.append(star(x0 + xm, Y_MET, r, BURST, LOSS, op=op, sw=LW_STAR))

    # --- emitter
    g.append(star(x0 + X_EMIT, Y_EMIT, 8.5, EMIT, EMIT_EDGE, n=9, inner=0.42, sw=LW_STAR))
    g.append('<circle cx="%s" cy="%s" r="2.4" fill="#fff8e6"/>'
             % (f(x0 + X_EMIT), f(Y_EMIT)))

    # --- titles
    g.append(txt(x0 + PW / 2, 28, title, 21, INK, "middle", "700"))
    g.append(txt(x0 + PW / 2, 52, sub, 16, MUTED, "middle"))
    g.append('</g>')
    return "\n".join(g), sum(esc)

# ---------------------------------------------------------------- assemble
A('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
  % (W, H, W, H))
A('<defs><linearGradient id="fade" x1="0" y1="0" x2="1" y2="1">'
  '<stop offset="0" stop-color="%s" stop-opacity="0.9"/>'
  '<stop offset="1" stop-color="%s" stop-opacity="0"/></linearGradient></defs>' % (GREEN, GREEN))
A('<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H))
A('<title>Repeated outcoupling attempts in a conventional OLED versus one designed '
  'with the proposed design rule</title>')

p1, e1 = panel(PANEL_X[0], .30, .72, .90, "Conventional OLED + outcoupling structure",
               "large absorption per round trip — the beam dies out within a few passes")
p2, e2 = panel(PANEL_X[1], .30, .98, .995, "Designed by the proposed design rule",
               "absorption suppressed — intensity survives many round trips")
A(p1); A(p2)

# ---------------------------------------------------------------- layer keys
for y_lab, y_tip, s in ((150, 152, "Outcoupling structure"),
                        (222, 218, "Glass substrate"),
                        (263, 269, "TCO anode"),
                        (292, 286, "Organic layers"),
                        (318, 315, "Metal cathode")):
    A(txt(181, y_lab + 5, s, 14, MUTED, "end"))
    A('<path d="M 187,%s L 195,%s" stroke="%s" stroke-width="0.9" fill="none"/>'
      % (LEADER, f(y_lab), f(y_tip)))

# ---------------------------------------------------------------- legend
LG = 412.0
A(ray((642, LG), (706, LG), 6.2))
A(txt(720, LG + 5, "Light ray  (line width ∝ optical power)", 16, MUTED))
A(star(1070, LG - 2, 8.5, EMIT, EMIT_EDGE, n=9, inner=0.42, sw=LW_STAR))
A(txt(1088, LG + 5, "Exciton emission", 15, MUTED))
A(wave(1300, LG - 2, 56, 3.4, 2.5, LOSS, 2.2, 0.95))
A(txt(1376, LG + 5, "Absorption loss  (ohmic at the metal, TCO)", 15, MUTED))
A('</svg>')

dest = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outcoupling_roundtrip.svg")
with open(dest, "w") as fh:
    fh.write("\n".join(out))
print("wrote %s   extracted: conventional %.0f%%  vs  design rule %.0f%%"
      % (dest, 100 * e1, 100 * e2))
