# -*- coding: utf-8 -*-
"""Generate an isometric exploded-stack figure of a tandem blue PEP-OLED as SVG."""

# ---------------------------------------------------------------- projection
OX = 22.0                 # x of the front-left corner of every top face
EX, EY = 330.0, 36.0      # "east" edge vector  (left -> right, sloping down)
DX, DY = 62.0, -72.0      # "depth" vector      (front -> back, up & right)
SLANT = 6.2               # degrees, atan(EY/EX) -- labels follow this

W, H = 440, 520

def f(v):
    s = "%.2f" % v
    return s.rstrip("0").rstrip(".") if "." in s else s

def pts(seq):
    return " ".join("%s,%s" % (f(x), f(y)) for x, y in seq)

def top_face(oy):
    """Parallelogram of the top face whose front-left corner is (OX, oy)."""
    return [(OX, oy), (OX + EX, oy + EY),
            (OX + EX + DX, oy + EY + DY), (OX + DX, oy + DY)]

def front_band(oy, t1, t2):
    """Slanted band on the front face, t1..t2 below the top face."""
    return [(OX, oy + t1), (OX + EX, oy + EY + t1),
            (OX + EX, oy + EY + t2), (OX, oy + t2)]

def side_band(oy, t1, t2):
    """Matching band on the right-hand side face."""
    x, y = OX + EX, oy + EY
    return [(x, y + t1), (x + DX, y + DY + t1),
            (x + DX, y + DY + t2), (x, y + t2)]

out = []
A = out.append

# ------------------------------------------------------------------- layers
# (key, oy, thickness, top-face fill, [(t1, t2, fill), ...] front strata)
GLASSY = "url(#gGlass)"

layers = [
    dict(key="anode", oy=406, h=42, top="url(#gAnodeTop)",
         strata=[(0, 9,  "url(#gItoA)"),
                 (9, 25, "url(#gAg)"),
                 (25, 28, "url(#gGold)"),
                 (28, 42, GLASSY)]),
    dict(key="htl", oy=380, h=26, top="url(#gHtlTop)",
         strata=[(0, 2.5, "url(#gVioEdge)"), (2.5, 26, "url(#gVioFront)")]),
    dict(key="eml_b", oy=294, h=30, top="url(#gEmlTop)",
         strata=[(0, 30, "url(#gEmlFront)")]),
    dict(key="cgl", oy=262, h=32, top="url(#gCglTop)",
         strata=[(0, 32, "url(#gCglFront)")]),
    dict(key="eml_t", oy=232, h=30, top="url(#gEmlTop)",
         strata=[(0, 30, "url(#gEmlFront)")]),
    dict(key="etl", oy=150, h=26, top="url(#gHtlTop)",
         strata=[(0, 2.5, "url(#gVioEdge)"), (2.5, 26, "url(#gVioFront)")]),
    dict(key="cath", oy=112, h=38, top="url(#gCathTop)",
         strata=[(0, 12, "url(#gCathA)"),
                 (12, 25, "url(#gAl)"),
                 (25, 38, "url(#gCathC)")]),
]
L = {d["key"]: d for d in layers}

def slab(d):
    oy, h = d["oy"], d["h"]
    s = ['<g>']
    s.append('<polygon points="%s" fill="%s"/>' % (pts(top_face(oy)), d["top"]))
    for t1, t2, fill in d["strata"]:
        s.append('<polygon points="%s" fill="%s"/>' % (pts(front_band(oy, t1, t2)), fill))
    for t1, t2, fill in d["strata"]:
        s.append('<polygon points="%s" fill="%s"/>' % (pts(side_band(oy, t1, t2)), fill))
    # shade the right-hand face as one piece
    s.append('<polygon points="%s" fill="#0b1020" opacity="0.20"/>' % pts(side_band(oy, 0, h)))
    # crisp edges
    s.append('<polygon points="%s" fill="none" stroke="#ffffff" stroke-opacity="0.35" stroke-width="0.8"/>'
             % pts(top_face(oy)))
    s.append('<path d="M %s L %s L %s" fill="none" stroke="#0b1020" stroke-opacity="0.16" stroke-width="0.7"/>'
             % (f(OX) + "," + f(oy + h),
                f(OX + EX) + "," + f(oy + EY + h),
                f(OX + EX + DX) + "," + f(oy + EY + DY + h)))
    s.append('</g>')
    return "\n".join(s)

def label(d, text, fill):
    """Right-aligned label riding on the front face, tilted to match it."""
    x = 344.0
    y = d["oy"] + (x - OX) * EY / EX + d["h"] / 2.0 + 4.6
    return ('<text x="%s" y="%s" transform="rotate(%s %s %s)" text-anchor="end" '
            'fill="%s" font-size="13.5" font-family="%s" letter-spacing="0.2">%s</text>'
            % (f(x), f(y), SLANT, f(x), f(y), fill, FONT, text))

def glow(cy_layer_oy, dy=0.0, soft=(152, 42), core=(96, 22), core_op=0.95):
    """Emissive ellipse centred on a layer face, optionally pushed downwards."""
    cx = OX + EX / 2.0 + DX / 2.0
    cy = cy_layer_oy + EY / 2.0 + DY / 2.0 + dy
    g = ['<g transform="rotate(%s %s %s)">' % (SLANT, f(cx), f(cy))]
    g.append('<ellipse cx="%s" cy="%s" rx="%s" ry="%s" fill="url(#gGlowSoft)" filter="url(#blurSoft)"/>'
             % (f(cx), f(cy), f(soft[0]), f(soft[1])))
    g.append('<ellipse cx="%s" cy="%s" rx="%s" ry="%s" fill="url(#gGlowCore)" opacity="%s" '
             'filter="url(#blurCore)"/>'
             % (f(cx), f(cy), f(core[0]), f(core[1]), core_op))
    g.append('</g>')
    return "\n".join(g)

FONT = "Helvetica Neue, Helvetica, Arial, Liberation Sans, sans-serif"

# --------------------------------------------------------------------- defs
defs = """
<defs>
  <linearGradient id="gAnodeTop" x1="0" y1="0" x2="1" y2="0.35">
    <stop offset="0"    stop-color="#ffffff"/>
    <stop offset="0.45" stop-color="#eef1f4"/>
    <stop offset="1"    stop-color="#f9fbfc"/>
  </linearGradient>
  <linearGradient id="gItoA" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"   stop-color="#e9edf0"/>
    <stop offset="0.5" stop-color="#f4f6f8"/>
    <stop offset="1"   stop-color="#e2e7eb"/>
  </linearGradient>
  <linearGradient id="gAg" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"    stop-color="#a8afb7"/>
    <stop offset="0.18" stop-color="#eceef1"/>
    <stop offset="0.42" stop-color="#8d949c"/>
    <stop offset="0.66" stop-color="#d3d8dd"/>
    <stop offset="0.86" stop-color="#8f959d"/>
    <stop offset="1"    stop-color="#b6bcc3"/>
  </linearGradient>
  <linearGradient id="gGold" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"   stop-color="#c9a24f"/>
    <stop offset="0.4" stop-color="#e8cd86"/>
    <stop offset="1"   stop-color="#c19a45"/>
  </linearGradient>
  <linearGradient id="gGlass" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"   stop-color="#fdfdfe"/>
    <stop offset="0.6" stop-color="#f2f4f6"/>
    <stop offset="1"   stop-color="#e9ecef"/>
  </linearGradient>

  <linearGradient id="gCathTop" x1="0" y1="0" x2="1" y2="0.35">
    <stop offset="0"    stop-color="#ffffff"/>
    <stop offset="0.5"  stop-color="#e9edf1"/>
    <stop offset="1"    stop-color="#fafbfc"/>
  </linearGradient>
  <linearGradient id="gCathA" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"   stop-color="#f6f8f9"/>
    <stop offset="0.5" stop-color="#ffffff"/>
    <stop offset="1"   stop-color="#eef1f3"/>
  </linearGradient>
  <linearGradient id="gAl" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"    stop-color="#d3d8de"/>
    <stop offset="0.22" stop-color="#eef1f4"/>
    <stop offset="0.5"  stop-color="#bfc5cd"/>
    <stop offset="0.78" stop-color="#e6e9ed"/>
    <stop offset="1"    stop-color="#c6ccd3"/>
  </linearGradient>
  <linearGradient id="gCathC" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"   stop-color="#fbfcfd"/>
    <stop offset="0.6" stop-color="#f0f2f5"/>
    <stop offset="1"   stop-color="#e7eaee"/>
  </linearGradient>

  <linearGradient id="gHtlTop" x1="0" y1="0" x2="1" y2="0.35">
    <stop offset="0"    stop-color="#8b46e6"/>
    <stop offset="0.45" stop-color="#a865f2"/>
    <stop offset="1"    stop-color="#7c3aed"/>
  </linearGradient>
  <linearGradient id="gVioEdge" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0"   stop-color="#c99cff"/>
    <stop offset="0.5" stop-color="#e0c4ff"/>
    <stop offset="1"   stop-color="#bb8cf8"/>
  </linearGradient>
  <linearGradient id="gVioFront" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0"   stop-color="#9333ea"/>
    <stop offset="0.55" stop-color="#7526d6"/>
    <stop offset="1"   stop-color="#5b1ab0"/>
  </linearGradient>

  <linearGradient id="gEmlTop" x1="0" y1="0" x2="1" y2="0.35">
    <stop offset="0"    stop-color="#2f5ad8"/>
    <stop offset="0.45" stop-color="#5386f5"/>
    <stop offset="1"    stop-color="#3462e0"/>
  </linearGradient>
  <linearGradient id="gEmlFront" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0"   stop-color="#2a4fc4"/>
    <stop offset="1"   stop-color="#17308c"/>
  </linearGradient>

  <linearGradient id="gCglTop" x1="0" y1="0" x2="1" y2="0.35">
    <stop offset="0"    stop-color="#aadeec"/>
    <stop offset="0.45" stop-color="#c8ecf5"/>
    <stop offset="1"    stop-color="#a3d9e9"/>
  </linearGradient>
  <linearGradient id="gCglFront" x1="0" y1="0" x2="0" y2="1">
    <stop offset="0"   stop-color="#7cc4d8"/>
    <stop offset="1"   stop-color="#529fba"/>
  </linearGradient>

  <radialGradient id="gGlowSoft" cx="0.5" cy="0.5" r="0.5">
    <stop offset="0"    stop-color="#ffffff" stop-opacity="0.95"/>
    <stop offset="0.32" stop-color="#8fe8ff" stop-opacity="0.80"/>
    <stop offset="0.68" stop-color="#3fb0ff" stop-opacity="0.34"/>
    <stop offset="1"    stop-color="#2a7fe0" stop-opacity="0"/>
  </radialGradient>
  <radialGradient id="gGlowCore" cx="0.5" cy="0.5" r="0.5">
    <stop offset="0"    stop-color="#ffffff" stop-opacity="1"/>
    <stop offset="0.5"  stop-color="#e4faff" stop-opacity="0.92"/>
    <stop offset="1"    stop-color="#9ceaff" stop-opacity="0"/>
  </radialGradient>
  <radialGradient id="gShadow" cx="0.5" cy="0.5" r="0.5">
    <stop offset="0"   stop-color="#2b3340" stop-opacity="0.26"/>
    <stop offset="0.6" stop-color="#2b3340" stop-opacity="0.10"/>
    <stop offset="1"   stop-color="#2b3340" stop-opacity="0"/>
  </radialGradient>

  <filter id="blurSoft" x="-50%" y="-120%" width="200%" height="340%">
    <feGaussianBlur stdDeviation="6"/>
  </filter>
  <filter id="blurCore" x="-50%" y="-120%" width="200%" height="340%">
    <feGaussianBlur stdDeviation="3"/>
  </filter>
  <filter id="blurShadow" x="-50%" y="-150%" width="200%" height="400%">
    <feGaussianBlur stdDeviation="5"/>
  </filter>
</defs>
"""

# ------------------------------------------------------------------ assemble
A('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
  'width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H))
A(defs)
A('<rect width="%d" height="%d" fill="#ffffff"/>' % (W, H))

A('<title>Tandem blue PEP-OLED device stack</title>')
A('<text x="220" y="26" text-anchor="middle" fill="#3c4149" font-size="15" '
  'font-family="%s">Tandem blue PEP-OLED</text>' % FONT)

# ground shadow
A('<g id="ground-shadow" transform="rotate(%s 222 494)">'
  '<ellipse cx="222" cy="494" rx="176" ry="20" fill="url(#gShadow)" filter="url(#blurShadow)"/></g>' % SLANT)

A('<g id="bottom-unit">')
A(slab(L["anode"]))
A(slab(L["htl"]))
A('</g>')

A('<g id="glow-lower">')
A(glow(L["htl"]["oy"], dy=-2, soft=(150, 40), core=(96, 21), core_op=0.95))
A('</g>')

A('<g id="middle-unit">')
A(slab(L["eml_b"]))
A(slab(L["cgl"]))
A(slab(L["eml_t"]))
A('</g>')

A('<g id="glow-upper">')
A(glow(L["etl"]["oy"] + L["etl"]["h"], dy=34, soft=(150, 40), core=(92, 20), core_op=0.9))
A('</g>')

A('<g id="top-unit">')
A(slab(L["etl"]))
A(slab(L["cath"]))
A('</g>')

A('<g id="labels">')
A(label(L["cath"],  "Ag/Al/Liq cathode", "#2b3038"))
A(label(L["etl"],   "Polaritonic ETL",   "#ffffff"))
A(label(L["eml_t"], "Blue EML",          "#ffffff"))
A(label(L["cgl"],   "CGL",               "#ffffff"))
A(label(L["eml_b"], "Blue EML",          "#ffffff"))
A(label(L["htl"],   "Polaritonic HTL",   "#ffffff"))
A(label(L["anode"], "ITO/Ag/ITO anode",  "#2b3038"))
A('</g>')

A('</svg>')

import os
dest = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "tandem_blue_pep_oled.svg")
with open(dest, "w") as fh:
    fh.write("\n".join(out))
print("wrote " + dest)
