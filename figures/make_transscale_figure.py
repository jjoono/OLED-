# -*- coding: utf-8 -*-
"""Trans-scale optical simulation schematic: cut-away isometric OLED stack."""
import math, os

# ------------------------------------------------------------------ projection
EX, EY = 400.0, 46.0        # "east"  : left -> right, sloping down
DX, DY = 104.0, -118.0      # "depth" : front -> back, up & right
OX0 = 70.0                  # front-left corner x of the s = 1 footprint
U0, V0 = 0.42, 0.55         # notch corner, in s = 1 footprint coordinates
W, H = 1020, 565

ELEN = math.hypot(EX, EY)
EUX, EUY = EX / ELEN, EY / ELEN
SLANT = math.degrees(math.atan2(EY, EX))
FONT = "Helvetica Neue, Helvetica, Arial, Liberation Sans, sans-serif"

def f(v):
    s = "%.2f" % v
    return s.rstrip("0").rstrip(".") if "." in s else s

def pts(seq):
    return " ".join("%s,%s" % (f(x), f(y)) for x, y in seq)

# ------------------------------------------------------------------- geometry
def frame(T, s):
    """Footprint of a layer: top plane T, footprint scaled by s about the centre."""
    return dict(O=(OX0 + (1 - s) * (EX + DX) / 2.0, T + (1 - s) * (EY + DY) / 2.0),
                E=(s * EX, s * EY), D=(s * DX, s * DY),
                u0=0.5 + (U0 - 0.5) / s, v0=0.5 + (V0 - 0.5) / s, T=T, s=s)

def P(fr, u, v, t=0.0):
    O, E, D = fr["O"], fr["E"], fr["D"]
    return (O[0] + u * E[0] + v * D[0], O[1] + u * E[1] + v * D[1] + t)

def top_hex(fr):
    u0, v0 = fr["u0"], fr["v0"]
    return [P(fr, 0, 0), P(fr, u0, 0), P(fr, u0, v0),
            P(fr, 1, v0), P(fr, 1, 1), P(fr, 0, 1)]

def wall(fr, a, b, h):
    return [P(fr, a[0], a[1]), P(fr, b[0], b[1]),
            P(fr, b[0], b[1], h), P(fr, a[0], a[1], h)]

def walls(fr, h):
    """front-outer, cut-east, cut-front, right-outer  (drawn back -> front)."""
    u0, v0 = fr["u0"], fr["v0"]
    return dict(front=wall(fr, (0, 0), (u0, 0), h),
                cut_e=wall(fr, (u0, 0), (u0, v0), h),
                cut_f=wall(fr, (u0, v0), (1, v0), h),
                right=wall(fr, (1, v0), (1, 1), h))

# ------------------------------------------------------------------- drawing
out = []
A = out.append

def poly(p, fill, extra=""):
    return '<polygon points="%s" fill="%s"%s/>' % (pts(p), fill, extra)

EDGE = ' stroke="#1d232b" stroke-opacity="0.16" stroke-width="0.7" stroke-linejoin="round"'

def slab(fr, h, key, lens=False):
    w = walls(fr, h)
    s = ['<g id="layer-%s">' % key]
    s.append(poly(top_hex(fr), "url(#g%sTop)" % key, EDGE))
    for name, grad in (("right", "Side"), ("cut_e", "Side"),
                       ("cut_f", "Front"), ("front", "Front")):
        s.append(poly(w[name], "url(#g%s%s)" % (key, grad), EDGE))
    # contact shading where the layer above sits on this one
    for name in ("cut_f", "front"):
        q = w[name]
        top = [q[0], q[1],
               (q[1][0], q[1][1] + 9), (q[0][0], q[0][1] + 9)]
        s.append(poly(top, "url(#gAO)"))
    if lens:                                   # micro-lens array on the underside
        for name, grad in (("front", "Front"), ("cut_f", "Front"),
                           ("cut_e", "Side"), ("right", "Side")):
            q = w[name]
            s.append(scallop_path(q[3], q[2], 11.0, "url(#g%s%s)" % (key, grad)))
    s.append('</g>')
    return "\n".join(s)

def scallop_path(p0, p1, bump, fill):
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    n = max(1, int(round(math.hypot(dx, dy) / 19.0)))
    k = bump * 1.34
    d = "M %s,%s" % (f(p0[0]), f(p0[1]))
    for i in range(n):
        a = (p0[0] + dx * i / n,       p0[1] + dy * i / n)
        b = (p0[0] + dx * (i + 1) / n, p0[1] + dy * (i + 1) / n)
        d += " C %s,%s %s,%s %s,%s" % (f(a[0]), f(a[1] + k), f(b[0]), f(b[1] + k),
                                       f(b[0]), f(b[1]))
    d += " Z"
    return '<path d="%s" fill="%s"%s/>' % (d, fill, EDGE)

def arrow(p0, p1, color, width=2.4, dash="8 6", head=12, op=1.0):
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    L = math.hypot(dx, dy) or 1.0
    ux, uy = dx / L, dy / L
    q = (p1[0] - ux * head, p1[1] - uy * head)
    px, py = -uy, ux
    hw = head * 0.42
    tri = [p1, (q[0] + px * hw, q[1] + py * hw), (q[0] - px * hw, q[1] - py * hw)]
    da = ' stroke-dasharray="%s"' % dash if dash else ""
    return ('<g opacity="%s"><line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" '
            'stroke-width="%s" stroke-linecap="round"%s/>'
            '<polygon points="%s" fill="%s"/></g>'
            % (op, f(p0[0]), f(p0[1]), f(q[0]), f(q[1]), color, f(width), da,
               pts(tri), color))

def arcpath(cx, cy, r, a0, a1, n=28):
    p = [(cx + r * math.cos(math.radians(a)), cy + r * math.sin(math.radians(a)))
         for a in [a0 + (a1 - a0) * i / n for i in range(n + 1)]]
    return "M " + " L ".join("%s,%s" % (f(x), f(y)) for x, y in p)

def brace(x, y0, y1, w=11.0, r=13.0, color="#8b9299"):
    ym = (y0 + y1) / 2.0
    d = ("M %s,%s Q %s,%s %s,%s L %s,%s Q %s,%s %s,%s "
         "Q %s,%s %s,%s L %s,%s Q %s,%s %s,%s"
         % (f(x), f(y0), f(x + w), f(y0), f(x + w), f(y0 + r),
            f(x + w), f(ym - r), f(x + w), f(ym), f(x + 2 * w), f(ym),
            f(x + w), f(ym), f(x + w), f(ym + r),
            f(x + w), f(y1 - r), f(x + w), f(y1), f(x), f(y1)))
    return ('<path d="%s" fill="none" stroke="%s" stroke-width="2.3" '
            'stroke-linecap="round" stroke-linejoin="round"/>' % (d, color))

def face_label(fr, h, lines, fill, size=13.0, pad=15.0, lead=17.0):
    O = fr["O"]
    x = O[0] + pad
    ymid = O[1] + pad * EY / EX + h / 2.0
    y0 = ymid - lead * (len(lines) - 1) / 2.0 + size * 0.35
    s = []
    for i, t in enumerate(lines):
        y = y0 + i * lead
        s.append('<text x="%s" y="%s" transform="rotate(%s %s %s)" fill="%s" '
                 'font-size="%s" font-family="%s">%s</text>'
                 % (f(x), f(y), f(SLANT), f(x), f(y), fill, f(size), FONT, t))
    return "\n".join(s)

# --------------------------------------------------------------- stack layout
LAYERS = [   # key, top plane, thickness, footprint scale
    ("Out",  278.0, 90.0, 1.14),
    ("Glass",234.0, 44.0, 1.14),
    ("Ito",  224.0, 10.0, 0.97),
    ("Oled", 152.0, 72.0, 1.00),
    ("Cath", 130.0, 22.0, 0.86),
]
FR = {k: frame(T, s) for k, T, h, s in LAYERS}
TH = {k: h for k, T, h, s in LAYERS}

# cut-plane local frame: a = along E from the notch corner, b = depth below T=130
NX = OX0 + U0 * EX + V0 * DX
NY = 130.0 + U0 * EY + V0 * DY
CUT = 'matrix(%s,%s,0,1,%s,%s)' % (f(EUX), f(EUY), f(NX), f(NY))
B = {k: (T - 130.0, T - 130.0 + h) for k, T, h, s in LAYERS}   # band depths

# ------------------------------------------------------------------ gradients
def lg(i, stops, x1=0, y1=0, x2=0, y2=1):
    s = ['<linearGradient id="%s" x1="%s" y1="%s" x2="%s" y2="%s">' % (i, x1, y1, x2, y2)]
    for off, col in stops:
        s.append('<stop offset="%s" stop-color="%s"/>' % (off, col))
    s.append('</linearGradient>')
    return "".join(s)

DEFS = ["<defs>"]
DEFS += [
    lg("gCathTop",  [(0, "#eef1f5"), (.35, "#c6ccd4"), (.62, "#e6eaef"), (1, "#cbd1d8")], 0, 0, 1, .35),
    lg("gCathFront",[(0, "#ced4db"), (1, "#a8b0ba")]),
    lg("gCathSide", [(0, "#9aa2ac"), (1, "#7f8892")]),
    lg("gOledTop",  [(0, "#ffd489"), (.45, "#fbb43f"), (1, "#f5a623")], 0, 0, 1, .35),
    lg("gOledFront",[(0, "#f2a91f"), (1, "#d5860c")]),
    lg("gOledSide", [(0, "#c47c09"), (1, "#a66407")]),
    lg("gItoTop",   [(0, "#7cbcff"), (1, "#3d8ef7")], 0, 0, 1, .35),
    lg("gItoFront", [(0, "#3d8ef7"), (1, "#2a6ed0")]),
    lg("gItoSide",  [(0, "#2560b8"), (1, "#1c4d96")]),
    lg("gGlassTop", [(0, "#cdeefb"), (.34, "#ffffff"), (.52, "#bfe9f8"), (1, "#d6f1fb")], 0, 0, 1, .35),
    lg("gGlassFront",[(0, "#aee2f4"), (1, "#84cde8")]),
    lg("gGlassSide",[(0, "#7cc3df"), (1, "#61accc")]),
    lg("gOutTop",   [(0, "#f7f9fa"), (.5, "#e7ebee"), (1, "#f2f5f6")], 0, 0, 1, .35),
    lg("gOutFront", [(0, "#edf0f2"), (1, "#d5dbe0")]),
    lg("gOutSide",  [(0, "#c8cfd6"), (1, "#b1b9c2")]),
    '<linearGradient id="gAO" x1="0" y1="0" x2="0" y2="1">'
    '<stop offset="0" stop-color="#0d1218" stop-opacity="0.17"/>'
    '<stop offset="1" stop-color="#0d1218" stop-opacity="0"/></linearGradient>',
    '<radialGradient id="gEmit" cx="0.5" cy="0.5" r="0.5">'
    '<stop offset="0" stop-color="#ffffff" stop-opacity="0.95"/>'
    '<stop offset="0.35" stop-color="#bfe6ff" stop-opacity="0.55"/>'
    '<stop offset="1" stop-color="#7fc4ff" stop-opacity="0"/></radialGradient>',
    '<radialGradient id="gSpill" cx="0.5" cy="0.35" r="0.65">'
    '<stop offset="0" stop-color="#cfe9ff" stop-opacity="0.62"/>'
    '<stop offset="0.45" stop-color="#bcdefb" stop-opacity="0.34"/>'
    '<stop offset="0.74" stop-color="#a9d4f7" stop-opacity="0.13"/>'
    '<stop offset="1" stop-color="#9cc8f0" stop-opacity="0"/></radialGradient>',
    '<radialGradient id="gBead" cx="0.34" cy="0.30" r="0.78">'
    '<stop offset="0" stop-color="#d5ecff"/><stop offset="0.45" stop-color="#5aa9f5"/>'
    '<stop offset="1" stop-color="#2a6fc4"/></radialGradient>',
]
# clip: the union of every layer's cut-away front face  +  per-regime bands
DEFS.append('<clipPath id="cutA">')
for k, T, h, s in LAYERS:
    DEFS.append('<polygon points="%s"/>' % pts(walls(FR[k], TH[k])["cut_f"]))
DEFS.append('</clipPath>')
for k, (b0, b1) in B.items():
    DEFS.append('<clipPath id="z%s"><rect x="-40" y="%s" width="460" height="%s"/></clipPath>'
                % (k, f(b0), f(b1 - b0)))
DEFS.append('<clipPath id="zWave"><rect x="-40" y="22" width="460" height="82"/></clipPath>')
DEFS.append("</defs>")

# ------------------------------------------------------------------- assemble
A('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">'
  % (W, H, W, H))
A("\n".join(DEFS))
A('<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H))
A('<title>Trans-scale optical simulation of an OLED with an outcoupling structure</title>')

# light spilling out of the bottom of the device
A('<ellipse cx="215" cy="424" rx="250" ry="58" fill="url(#gSpill)"/>')

for k, T, h, s in LAYERS:
    A(slab(FR[k], h, k, lens=(k == "Out")))

# ------------------------------------------------- physics on the cut plane
A('<g clip-path="url(#cutA)"><g transform="%s">' % CUT)

DIP = (118.0, 58.0)
AX = math.radians(125.0)
dxu, dyu = math.cos(AX), math.sin(AX)
pxu, pyu = -dyu, dxu

# spreading wavefronts inside the organic stack
A('<g clip-path="url(#zWave)" fill="none" stroke="#1f6fb2" stroke-linecap="round">')
for r, op in ((34, .46), (48, .34), (62, .23), (76, .14)):
    A('<path d="%s" stroke-width="1.6" stroke-opacity="%s"/>'
      % (arcpath(DIP[0], DIP[1], r, 8, 172), op))
A('</g>')

# dipole emission lobes (two circles tangent at the emitter) + dipole moment
A('<g clip-path="url(#zOled)">')
A('<ellipse cx="%s" cy="%s" rx="46" ry="40" fill="url(#gEmit)"/>' % (f(DIP[0]), f(DIP[1])))
for r in (20.0, 11.5):
    for sgn in (1, -1):
        A('<circle cx="%s" cy="%s" r="%s" fill="none" stroke="#3f4a57" '
          'stroke-width="1.5" stroke-opacity="0.85"/>'
          % (f(DIP[0] + sgn * r * pxu), f(DIP[1] + sgn * r * pyu), f(r)))
A(arrow((DIP[0] - 30 * dxu, DIP[1] - 30 * dyu),
        (DIP[0] + 34 * dxu, DIP[1] + 34 * dyu), "#14181e", 2.4, "", 11))
A('<circle cx="%s" cy="%s" r="4.6" fill="#ffffff"/>' % (f(DIP[0]), f(DIP[1])))
A('<circle cx="%s" cy="%s" r="2.9" fill="#14181e"/>' % (f(DIP[0]), f(DIP[1])))
A('</g>')

# rays crossing the glass
A('<g clip-path="url(#zGlass)">')
for p0, p1 in (((60, 106), (34, 145)), ((114, 106), (114, 146)), ((168, 106), (196, 145))):
    A(arrow(p0, p1, "#20262e", 2.3, "7 5", 11))
A('</g>')

# outcoupling layer: far-field lobe on the left, scattering beads on the right
A('<g clip-path="url(#zOut)">')
CX, CY, R = 64.0, 168.0, 41.0
A('<g fill="none" stroke="#9aa3ac" stroke-width="0.9" stroke-opacity="0.75">')
for rr in (R, R * 0.62, R * 0.3):
    A('<path d="%s"/>' % arcpath(CX, CY, rr, 0, 180))
for ang in (0, 45, 90, 135, 180):
    A('<line x1="%s" y1="%s" x2="%s" y2="%s"/>'
      % (f(CX), f(CY), f(CX + R * math.cos(math.radians(ang))),
         f(CY + R * math.sin(math.radians(ang)))))
A('</g>')
lobe = []
for i in range(97):
    th = math.radians(i * 180.0 / 96.0)
    rr = R * (0.30 + 0.70 * ((1 + math.cos(th - math.radians(90))) / 2.0) ** 1.25)
    lobe.append((CX + rr * math.cos(th), CY + rr * math.sin(th)))
A('<path d="M %s Z" fill="#e2342e" fill-opacity="0.10" stroke="#e2342e" '
  'stroke-width="2.2" stroke-linejoin="round"/>'
  % " L ".join("%s,%s" % (f(x), f(y)) for x, y in lobe))
A(arrow((28, 120), (CX - 2, CY - 4), "#20262e", 2.3, "7 5", 11))

BEADS = [(138, 170, 7.0), (174, 159, 6.0), (210, 172, 8.0), (152, 199, 6.5),
         (192, 203, 7.5), (230, 193, 6.0), (245, 165, 5.5), (170, 220, 6.0),
         (218, 222, 5.2), (256, 205, 6.4)]
A('<g stroke="#20262e" stroke-opacity="0.45" stroke-width="1.3" '
  'stroke-dasharray="5 4" stroke-linecap="round">')
for a, b in (((120, 150), (135, 165)), ((141, 175), (170, 216)),
             ((143, 173), (205, 168)), ((214, 177), (228, 190)),
             ((196, 208), (215, 218))):
    A('<line x1="%s" y1="%s" x2="%s" y2="%s"/>' % (f(a[0]), f(a[1]), f(b[0]), f(b[1])))
A('</g>')
for x, y, r in BEADS:
    A('<circle cx="%s" cy="%s" r="%s" fill="url(#gBead)"/>' % (f(x), f(y), f(r)))
A('</g>')
A('</g></g>')

# ------------------------------------------------------- outcoupled light
for p0, p1 in (((80, 390), (26, 470)), ((136, 397), (136, 482)), ((192, 403), (258, 470))):
    A(arrow(p0, p1, "#c62828", 4.0, "11 8", 17))

# ------------------------------------------------------------ layer captions
A(face_label(FR["Cath"],  TH["Cath"],  ["Metal cathode"],            "#232931", 15))
A(face_label(FR["Oled"],  TH["Oled"],  ["OLED"],                     "#3d2502", 19))
A(face_label(FR["Glass"], TH["Glass"], ["Glass"],                    "#10424d", 17))
A(face_label(FR["Out"],   TH["Out"],   ["Outcoupling", "structure"], "#2b3138", 16, 15.0, 20.0))

# ------------------------------------------------------------- regime panel
A(brace(640, 66, 282))
A(brace(640, 288, 442))
A('<rect x="690" y="66" width="290" height="376" rx="22" fill="none" '
  'stroke="#f2b705" stroke-width="3" stroke-dasharray="12 9"/>')
for y, t in ((182, "Wave-optics"),):
    A('<text x="835" y="%d" text-anchor="middle" fill="#1e2228" font-size="21" '
      'font-family="%s">%s</text>' % (y, FONT, t))
for y, t in ((357, "Geometrical"), (387, "(Ray) optics")):
    A('<text x="835" y="%d" text-anchor="middle" fill="#1e2228" font-size="21" '
      'font-family="%s">%s</text>' % (y, FONT, t))
for y, t in ((490, "Trans-scale"), (518, "optical simulation")):
    A('<text x="835" y="%d" text-anchor="middle" fill="#12161b" font-size="21" '
      'font-weight="700" font-family="%s">%s</text>' % (y, FONT, t))

A('</svg>')

dest = os.path.join(os.path.dirname(os.path.abspath(__file__)), "trans_scale_simulation.svg")
with open(dest, "w") as fh:
    fh.write("\n".join(out))
print("wrote " + dest)
