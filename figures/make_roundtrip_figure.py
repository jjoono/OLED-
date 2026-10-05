# -*- coding: utf-8 -*-
"""Fig.1(a): light recycling in an OLED with an external outcoupling structure.
Top    - conventional device: strong ohmic / TCO absorption, beam dies out.
Bottom - low-loss device: absorption suppressed, beam survives many round trips.

One figure.  Each row pairs the 3D render of a device with its cross-section: the
render shows how far the light actually spreads, the cross-section shows why.  The
layer names sit on the cross-section only -- carrying them twice is what made the
3D panel look like it was adding nothing.
"""
import math, os

W, H = 1620, 790
FONT = "Helvetica Neue, Helvetica, Arial, Liberation Sans, sans-serif"

# ---------------------------------------------------------------- stack (y)
Y_CREST, Y_TOP = 145.0, 175.0      # lens apex / lens base = glass top surface
Y_TCO, Y_ORG, Y_MET, Y_BOT = 261.0, 277.0, 295.0, 335.0
Y_EMIT = 286.0                     # emitting plane, middle of the organic stack
LENS_R, PITCH = 30.0, 60.0
PW = 900.0                         # panel width = 15 lenses
RENDER_X, RENDER_Y, RENDER_W = 55.0, 35.0, 300.0   # 3D render slot, square, per panel
PANEL_X = 640.0                    # both panels share a left edge; they stack
PANEL_DY = (0.0, 365.0)            # and are offset vertically instead.  A panel is
                                   # 335 tall: the escaping rays reach ~85 above the lens
                                   # crests, so the title has to sit clear of them

# Per panel, because the two devices differ exactly in these two layers.
TITLES = (("Conventional metal reflector",
           "large absorption per round trip \u2014 the beam dies out within a few passes"),
          ("Low loss reflector",
           "absorption suppressed \u2014 intensity survives many round trips"))
LAYER_Y = ((150, 152), (222, 218), (263, 269), (292, 286), (318, 315))
LAYER_KEYS = (("Outcoupling structure", "Glass substrate",
               "Transparent conductive oxide (TCO)", "Organic layers",
               "Conventional metal electrode"),
              ("Outcoupling structure", "Glass substrate",
               "Low loss TCO", "Organic layers",
               "Low loss metal electrode"))

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


def render_crop(tag):
    """Crop a device render to its own content and cache it beside the original.

    The renders are 1400x1400 with the device in the middle; dropped into a small
    slot whole, most of the slot is empty.  The box is squared off so the slot's
    aspect is preserved, and the file is rewritten only when it is missing or older
    than its source.
    """
    import os
    from PIL import Image
    here = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(here, "emitting_area_render_%s.png" % tag)
    dst = os.path.join(here, "emitting_area_render_%s_crop.png" % tag)
    if not os.path.exists(src):
        return None
    if os.path.exists(dst) and os.path.getmtime(dst) >= os.path.getmtime(src):
        return dst
    im = Image.open(src).convert("RGBA")
    # threshold first: the glow fades to alpha 1-2 across almost the whole
    # frame, so a plain getbbox() on alpha crops nothing
    bb = im.getchannel("A").point(lambda v: 255 if v > 10 else 0).getbbox()
    if bb:
        x0, y0, x1, y1 = bb
        side = max(x1 - x0, y1 - y0)
        cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        pad = side * 0.04
        side += 2 * pad
        im = im.crop((int(cx - side / 2), int(cy - side / 2),
                      int(cx + side / 2), int(cy + side / 2)))
    # 900 px is ~2x what the slot needs even at a 2400 px export, and it keeps the
    # base64 copy embedded in the SVG down to a few hundred kB
    if im.width > 900:
        im = im.resize((900, max(1, round(900 * im.height / im.width))), Image.LANCZOS)
    im.save(dst)
    return dst


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
A('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
  'width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H))
A('<defs><linearGradient id="fade" x1="0" y1="0" x2="1" y2="1">'
  '<stop offset="0" stop-color="%s" stop-opacity="0.9"/>'
  '<stop offset="1" stop-color="%s" stop-opacity="0"/></linearGradient></defs>' % (GREEN, GREEN))
A('<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H))
A('<title>Repeated outcoupling attempts in a conventional OLED versus one designed '
  'with the proposed design rule</title>')

esc = []
for i, (dy, (title, sub), keys, loss) in enumerate(
        zip(PANEL_DY, TITLES, LAYER_KEYS, ((.72, .90), (.98, .995)))):
    body, e = panel(PANEL_X, .30, loss[0], loss[1], title, sub)
    esc.append(e)
    A('<g transform="translate(0,%s)">' % f(dy))
    crop = render_crop(("conventional", "designrule")[i])
    if crop:
        # embedded, not referenced: a relative href breaks the moment the .svg is
        # moved, and most renderers refuse to load local files from an SVG at all
        import base64
        with open(crop, "rb") as fh:
            b64 = base64.b64encode(fh.read()).decode("ascii")
        A('<image xlink:href="data:image/png;base64,%s" x="%s" y="%s" '
          'width="%s" height="%s"/>'
          % (b64, f(RENDER_X), f(RENDER_Y), f(RENDER_W), f(RENDER_W)))
    A(body)
    for (y_lab, y_tip), s_lab in zip(LAYER_Y, keys):
        A(txt(PANEL_X - 19, y_lab + 5, s_lab, 14, MUTED, "end"))
        A('<path d="M %s,%s L %s,%s" stroke="%s" stroke-width="0.9" fill="none"/>'
          % (f(PANEL_X - 13), f(y_lab), f(PANEL_X - 5), f(y_tip), LEADER))
    A('</g>')
e1, e2 = esc

# ---------------------------------------------------------------- legend
LG = 745.0
A(ray((280, LG), (344, LG), 6.2))
A(txt(358, LG + 5, "Light ray  (line width \u221d optical power)", 16, MUTED))
A(star(720, LG - 2, 8.5, EMIT, EMIT_EDGE, n=9, inner=0.42, sw=LW_STAR))
A(txt(738, LG + 5, "Exciton emission", 15, MUTED))
A(wave(920, LG - 2, 56, 3.4, 2.5, LOSS, 2.2, 0.95))
A(txt(996, LG + 5, "Absorption loss  (ohmic at the metal, TCO)", 15, MUTED))
A('</svg>')

dest = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outcoupling_roundtrip.svg")
with open(dest, "w") as fh:
    fh.write("\n".join(out))
print("wrote %s   extracted: conventional %.0f%%  vs  design rule %.0f%%"
      % (dest, 100 * e1, 100 * e2))
