# -*- coding: utf-8 -*-
"""Two OLED layer stacks as an isometric PowerPoint figure, in the style of the
exploded-block device schematics Nature-family papers use.

Everything is a native PowerPoint shape: each layer is three freeform faces (top,
front-left, front-right), each label is a real text box, each emission arrow is a
straight arrow connector.  Nothing is grouped and nothing is a picture, so a layer
can be recoloured, renamed, resized or deleted in PowerPoint directly.

Layer order in STACKS is bottom -> top, the way the device is built.
"""
import math
import os

from pptx import Presentation
from pptx.util import Emu, Pt, Inches
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR, MSO_AUTO_SIZE
from pptx.oxml.ns import qn

# ------------------------------------------------------------------ canvas
W, H = 1400.0, 790.0                       # figure units
SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)
SCALE = SLIDE_W / W                        # EMU per figure unit
OX = 0
OY = int((SLIDE_H - H * SCALE) / 2.0)

def X(u): return int(round(OX + u * SCALE))
def Y(u): return int(round(OY + u * SCALE))
def D(u): return int(round(u * SCALE))
def PTS(u): return Pt(u * SCALE / 12700.0)

C = lambda h: RGBColor.from_string(h.lstrip("#").upper())

# ------------------------------------------------------------------ isometric
# screen x = 0.866 * (px - py)          px runs along the block's width
# screen y = 0.500 * (px + py) - pz     py along its depth, pz up
EX, EY = 0.866, 0.500
BW, BD = 150.0, 150.0                      # block width and depth: equal, so each layer
                                           # is a square film seen corner-on rather than a
                                           # rectangle.  BW + BD is unchanged, so the figure
                                           # keeps its screen footprint and the layout holds.
LH = 33.0                                  # default layer thickness on screen
INK, MUTED = C("#1b2026"), C("#5b646e")
RED = C("#d62828")                         # the source figures pick layers out in red
GREEN = C("#2e9e4f")
LEADER = C("#98a2ac")
LW = Emu(9525)                             # 0.75 pt outlines

def shade(hexcol, f):
    """Same hue, scaled toward black (f < 1) -- the three faces of one layer."""
    r, g, b = (int(hexcol[i:i + 2], 16) for i in (1, 3, 5))
    return RGBColor(*(min(255, int(round(v * f))) for v in (r, g, b)))

# ------------------------------------------------------------------ the stacks
# (label, fill, thickness, label colour, micro-lens array on its underside?)
GLASS, TCO = "#cfe4f5", "#7fd3e8"
HIL, HTL, EML = "#f6dfae", "#f0b860", "#b23a3a"
ETL, EIL = "#7f9ada", "#4560a8"
AG, ALU = "#c9ced3", "#b9bfc5"
NANO, PARY, REFL = "#d9e6e2", "#e2dcef", "#9aa1a8"
DBR_HI, DBR_LO = "#a02828", "#2b2f33"

STACK_A = [
    ("Glass + Microlens array film", GLASS, 46.0, RED, True),
    ("ITO [150]",                    TCO,   33.0, INK, False),
    ("HATCN/TAPC [5/55] 3 pairs",    HIL,   33.0, INK, False),
    ("TCTA [10]",                    HTL,   28.0, INK, False),
    (u"TCTA:B3PyMPM:Irppy\u2082acac [25]", EML, 36.0, INK, False),
    ("B3PyMPM [50]",                 ETL,   33.0, INK, False),
    (u"B3PyMPM:Cs\u2082CO\u2083 [5]", EIL, 26.0, INK, False),
    ("Ag [100]",                     AG,    30.0, INK, False),
    ("Nanolaminate",                 NANO,  33.0, INK, False),
    ("Parylene-C [3000]",            PARY,  40.0, RED, False),
    ("MoOx [5] + Ag reflector [100]", REFL, 34.0, RED, False),
]

STACK_B = [
    ("MLA substrate (n = 1.77)",     GLASS, 46.0, INK, True),
    ("IZO",                          TCO,   33.0, INK, False),
    ("TAPC/HAT-CN/TAPC/HAT-CN",      HIL,   33.0, INK, False),
    ("TCTA",                         HTL,   28.0, INK, False),
    (u"TCTA:B3PyMPM:Ir(dmppy-ph)\u2082tmd", EML, 36.0, INK, False),
    ("B3PyMPM",                      ETL,   33.0, INK, False),
    ("Al/Liq",                       ALU,   30.0, INK, False),
    ("DBR",                          None,  54.0, INK, False),   # alternating pairs
]

# A layer name may carry a second line: drawn smaller and in red under the label,
# for a layer that has more than one possible realisation.
STACK_C = [
    ("Glass (n = 1.77)",                                GLASS, 46.0, INK, False),
    ("Transparent electrode\nITO 50 nm   or   Ag 10 nm", TCO,  34.0, INK, False),
    ("HTL (200 nm)",                                    HTL,   40.0, INK, False),
    ("EML (20 nm)",                                     EML,   24.0, INK, False),
    ("ETL (200 nm)",                                    ETL,   40.0, INK, False),
    ("Ag (100 nm)",                                     AG,    30.0, INK, False),
]

DBR_PAIRS = 4                              # high/low index pairs drawn in the DBR block
MLA_PITCH = 19.0                           # micro-lens pitch on the substrate's underside
MLA_SEGS = 9                               # line segments per lens arc
FW = 300.0                                 # front-view: width of the front face
FDX, FDY = 32.0, 20.0                      # and how far back the top/side faces step
FMIN = 32.0                                # smallest layer height that still fits a label
FMIN2 = 50.0                               # ... and one that fits a label plus its note

SH = None                                  # set per figure by figure()


# ------------------------------------------------------------------ helpers
def _alpha(parent, a):
    if a is None or a >= 0.999 or parent is None:
        return
    sf = parent.find(qn("a:solidFill"))
    if sf is None:
        return
    clr = sf.find(qn("a:srgbClr"))
    clr.append(clr.makeelement(qn("a:alpha"), {"val": str(int(a * 100000))}))


def poly(pts, fill, line, name, lw=LW, op=None):
    """One flat face, as a closed freeform in figure coordinates."""
    dev = [(X(x), Y(y)) for x, y in pts]
    b = SH.build_freeform(dev[0][0], dev[0][1], scale=1.0)
    b.add_line_segments(dev[1:], close=True)
    s = b.convert_to_shape()
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    _alpha(s._element.spPr, op)
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line
        s.line.width = lw
    s.shadow.inherit = False
    s.name = name
    return s


def text(x, y, body, size, colour, align=PP_ALIGN.LEFT, bold=False,
         box_w=420.0, name="label", anchor=MSO_ANCHOR.MIDDLE):
    """`body` is a string, or a list of (text, subscript?) pairs."""
    left = {PP_ALIGN.LEFT: x, PP_ALIGN.CENTER: x - box_w / 2.0,
            PP_ALIGN.RIGHT: x - box_w}[align]
    h = size * 2.6
    tb = SH.add_textbox(X(left), Y(y - h / 2.0), D(box_w), D(h))
    tf = tb.text_frame
    tf.word_wrap = True                    # wrap="none" would override `align`
    tf.auto_size = MSO_AUTO_SIZE.NONE
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    p.alignment = align
    for chunk, sub in ([(body, False)] if isinstance(body, str) else body):
        r = p.add_run()
        r.text = chunk
        r.font.size = PTS(size)
        r.font.bold = bold
        r.font.name = "Arial"
        r.font.color.rgb = colour
        if sub:                            # a real PowerPoint subscript, still editable
            r.font._rPr.set("baseline", "-25000")
    tb.name = name
    return tb


def arrow(p0, p1, colour, w=2.6, name="emission", head=True):
    cn = SH.add_connector(MSO_CONNECTOR.STRAIGHT, X(p0[0]), Y(p0[1]), X(p1[0]), Y(p1[1]))
    cn.line.color.rgb = colour
    cn.line.width = Emu(max(3175, D(w)))
    ln = cn.line._get_or_add_ln()
    if head:
        ln.append(ln.makeelement(qn("a:tailEnd"),
                                 {"type": "triangle", "w": "med", "len": "med"}))
    cn.shadow.inherit = False
    cn.name = name
    return cn


def lens_profile(length):
    """A micro-lens array in cross-section: touching semicircles along `length`.

    Returns (t, drop) pairs, t running 0 -> length and drop hanging below the face.
    The pitch is adjusted to divide the edge exactly, so the array ends on a whole
    lens at both corners instead of a clipped one.
    """
    n = max(1, int(round(length / MLA_PITCH)))
    pitch = length / n
    r = pitch / 2.0
    out = []
    for i in range(n):
        c = (i + 0.5) * pitch
        for k in range(MLA_SEGS + 1):
            t = i * pitch + pitch * k / MLA_SEGS
            out.append((t, math.sqrt(max(0.0, r * r - (t - c) ** 2))))
    return out


def mla_depth(length):
    n = max(1, int(round(length / MLA_PITCH)))
    return length / n / 2.0


def slab(ox, oy, z0, z1, fill, tag, top=True, lens_bottom=False):
    """One layer: top face, front-left face (px = BW), front-right face (py = BD).

    Only the faces that can be seen are drawn -- an isometric box needs three, and
    the top one only on the layer that has nothing above it.

    `lens_bottom` scallops the two visible bottom edges into a micro-lens array,
    which is what the substrate of a bottom-emitting device actually looks like in
    section.  The bottom face itself is never visible from this angle, so the edges
    are the only place the array can show.
    """
    def P(px, py, pz):
        return (ox + EX * (px - py), oy + EY * (px + py) - pz)

    base = fill
    if top:
        poly([P(0, 0, z1), P(BW, 0, z1), P(BW, BD, z1), P(0, BD, z1)],
             shade(base, 1.00), C("#5d666f"), tag + " top")

    if lens_bottom:
        right = ([P(BW, 0, z1), P(BW, BD, z1)] +
                 [P(BW, BD - t, z0 - d) for t, d in lens_profile(BD)])
        front = ([P(BW, BD, z1), P(0, BD, z1)] +
                 [P(t, BD, z0 - d) for t, d in lens_profile(BW)])
    else:
        right = [P(BW, 0, z0), P(BW, 0, z1), P(BW, BD, z1), P(BW, BD, z0)]
        front = [P(BW, BD, z0), P(BW, BD, z1), P(0, BD, z1), P(0, BD, z0)]
    poly(right, shade(base, 0.86), C("#5d666f"), tag + " right")
    poly(front, shade(base, 0.70), C("#5d666f"), tag + " front")
    return P


# --------------------------------------------------------------- front view
def front_slab(x0, ybase, z0, z1, fill, tag, top=False, lens_bottom=False):
    """One layer seen head-on: front face, a sliver of the right side, and the top
    face on the layer that has nothing above it.

    The reference figures this imitates put the layer name inside the band, which
    only works head-on -- in the isometric view the front face is a sheared
    parallelogram and horizontal text sits on it badly.
    """
    F = lambda x, z: (x0 + x, ybase - z)
    B = lambda p: (p[0] + FDX, p[1] - FDY)
    if lens_bottom:
        bottom = [F(t, z0 - d) for t, d in lens_profile(FW)][::-1]
    else:
        bottom = [F(FW, z0), F(0, z0)]
    poly([F(0, z1), F(FW, z1)] + bottom, shade(fill, 1.00), C("#5d666f"), tag + " front")
    poly([F(FW, z0), B(F(FW, z0)), B(F(FW, z1)), F(FW, z1)],
         shade(fill, 0.78), C("#5d666f"), tag + " side")
    if top:
        poly([F(0, z1), F(FW, z1), B(F(FW, z1)), B(F(0, z1))],
             shade(fill, 1.10), C("#5d666f"), tag + " top")


def front_device(x0, ybase, layers, title, title_y):
    """`ybase` is the screen y of the stack's bottom edge."""
    z = 0.0
    for i, (name, fill, t, col, lens) in enumerate(layers):
        head, _, note = name.partition("\n")
        h = max(t, FMIN2 if note else FMIN)          # every band has to fit its label
        z0, z1, is_top = z, z + h, (i == len(layers) - 1)
        tag = "%s / %s" % (title, head)
        if fill is None:                             # the DBR: alternating pairs
            sub = h / (2.0 * DBR_PAIRS)
            for k in range(2 * DBR_PAIRS):
                front_slab(x0, ybase, z0 + k * sub, z0 + (k + 1) * sub,
                           DBR_HI if k % 2 == 0 else DBR_LO,
                           "%s pair %d" % (tag, k // 2 + 1),
                           top=is_top and k == 2 * DBR_PAIRS - 1)
            band = DBR_LO
        else:
            front_slab(x0, ybase, z0, z1, fill, tag, top=is_top, lens_bottom=lens)
            band = fill
        lum = sum(w * int(band[k:k + 2], 16) / 255.0
                  for w, k in ((0.2126, 1), (0.7152, 3), (0.0722, 5)))
        ink = C("#f2f5f8") if lum < 0.45 else col
        ym = ybase - (z0 + z1) / 2.0
        if note:
            text(x0 + FW / 2.0, ym - 11.0, _runs(head), 14.5, ink, PP_ALIGN.CENTER,
                 box_w=FW - 16.0, name=tag + " label")
            text(x0 + FW / 2.0, ym + 11.0, _runs(note), 11.5, C("#ffd7d7") if lum < 0.45 else RED,
                 PP_ALIGN.CENTER, box_w=FW - 16.0, name=tag + " label note")
        else:
            text(x0 + FW / 2.0, ym, _runs(head), 14.5, ink, PP_ALIGN.CENTER,
                 box_w=FW - 16.0, name=tag + " label")
        z = z1
    text(x0 + FW / 2.0 + FDX / 2.0, title_y, title, 19.0, INK, PP_ALIGN.CENTER,
         bold=True, box_w=520.0, name=title + " title")
    return z


def front_emission(x0, ybase, lenses):
    drop = mla_depth(FW) if lenses else 0.0
    cx = x0 + FW / 2.0
    for dx, dy in ((-86.0, 62.0), (0.0, 78.0), (86.0, 62.0)):
        arrow((cx + dx * 0.34, ybase + drop + 4), (cx + dx, ybase + drop + dy),
              GREEN, 2.6, "bottom emission")
    text(cx, ybase + drop + 112.0, "Bottom emission", 15.5, GREEN, PP_ALIGN.CENTER,
         box_w=300.0, name="emission label")


# ------------------------------------------------------------------ one device
def device(ox, oy, layers, title, label_x, label_w, title_y):
    """`oy` is the screen y of the stack's bottom-back corner (pz = 0)."""
    total = sum(t for _, _, t, _, _ in layers)
    z = 0.0
    right_x = ox + EX * BW                 # screen x of the block's right-hand corner
    for i, (name, fill, t, col, lens) in enumerate(layers):
        z0, z1, is_top = z, z + t, (i == len(layers) - 1)
        tag = "%s / %s" % (title, name.replace("\n", " / "))
        if fill is None:                   # the DBR: alternating quarter-wave pairs
            sub = t / (2.0 * DBR_PAIRS)
            for k in range(2 * DBR_PAIRS):
                slab(ox, oy, z0 + k * sub, z0 + (k + 1) * sub,
                     DBR_HI if k % 2 == 0 else DBR_LO,
                     "%s pair %d" % (tag, k // 2 + 1),
                     top=is_top and k == 2 * DBR_PAIRS - 1)
        else:
            slab(ox, oy, z0, z1, fill, tag, top=is_top, lens_bottom=lens)

        ym = oy + EY * BW - (z0 + z1) / 2.0        # the block's right-hand corner edge
        arrow((right_x + 6, ym), (label_x - 8, ym), LEADER, 1.1, tag + " leader", head=False)
        head, _, note = name.partition("\n")
        if note:
            text(label_x, ym - 12.0, _runs(head), 14.5, col, PP_ALIGN.LEFT,
                 box_w=label_w, name=tag + " label")
            text(label_x, ym + 12.0, _runs(note), 12.0, RED, PP_ALIGN.LEFT,
                 box_w=label_w, name=tag + " label note")
        else:
            text(label_x, ym, _runs(name), 14.5, col, PP_ALIGN.LEFT,
                 box_w=label_w, name=tag + " label")
        z = z1

    text(ox + EX * (BW - BD) / 2.0, title_y, title, 19.0, INK, PP_ALIGN.CENTER, bold=True,
         box_w=520.0, name=title + " title")
    return total


def _runs(name):
    """Split a layer name so digits after an element symbol become real subscripts."""
    out, buf = [], ""
    for ch in name:
        if ch in u"₂₃":
            if buf:
                out.append((buf, False)); buf = ""
            out.append(({u"₂": "2", u"₃": "3"}[ch], True))
        else:
            buf += ch
    if buf:
        out.append((buf, False))
    return out


def emission(ox, oy, lenses=True):
    """Green arrows leaving the substrate: both devices emit through the bottom.

    The tails sit on the two visible bottom edges, not at the centre of the bottom
    face -- that point is behind the block, and arrows starting there look like they
    come out of the middle of the stack.
    """
    dw, dd = (mla_depth(BW), mla_depth(BD)) if lenses else (0.0, 0.0)
    corner = (ox + EX * (BW - BD), oy + EY * (BW + BD) + max(dw, dd))
    left = (ox - EX * BD + 0.60 * EX * BW, oy + EY * BD + 0.60 * EY * BW + dw)
    right = (ox + EX * BW - 0.60 * EX * BD, oy + EY * BW + 0.60 * EY * BD + dd)
    for (tx, ty), (hx, hy) in ((left, (-58.0, 46.0)), (corner, (0.0, 68.0)), (right, (58.0, 46.0))):
        arrow((tx, ty + 3), (tx + hx, ty + hy), GREEN, 2.6, "bottom emission")
    text(corner[0], corner[1] + 100.0, "Bottom emission", 15.5, GREEN, PP_ALIGN.CENTER,
         box_w=300.0, name="emission label")


# ------------------------------------------------------------------ build
TITLE_GAP = 62.0                           # between a title and the top of its stack


def figure(panels, filename):
    """One slide, one file.  `panels` is (ox, layers, title, label_x, label_w).

    The baseline follows the tallest stack rather than being fixed: a six-layer
    stack under a baseline set for an eleven-layer one leaves the title marooned
    at the top of an empty slide.
    """
    global SH
    tallest = max(sum(t for _, _, t, _, _ in layers) for _, layers, _, _, _ in panels)
    lens = any(layers[0][4] for _, layers, _, _, _ in panels)
    below = EY * (BW + BD) + (mla_depth(BW) if lens else 0.0) + 122.0
    base_y = H / 2.0 + (tallest + TITLE_GAP + 25.0 - below) / 2.0
    title_y = base_y - tallest - TITLE_GAP

    prs = Presentation()
    prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H
    SH = prs.slides.add_slide(prs.slide_layouts[6]).shapes
    for ox, layers, title, label_x, label_w in panels:
        device(ox, base_y, layers, title, label_x, label_w, title_y)
        emission(ox, base_y, lenses=layers[0][4])
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    prs.save(out)
    print("saved %s  (%d shapes, 0 groups, 0 pictures)" % (out, len(SH)))


def front_figure(panels, filename, fw=None):
    """Front view: `panels` is (x0, layers, title).

    `fw` widens the stack for a figure that holds only one device, where the
    default width leaves it marooned in the middle of the slide.  Text in a
    PowerPoint shape does not scale when the shape is resized, so this has to be
    right at build time rather than left to the reader.
    """
    global SH, FW
    FW = fw or 300.0
    def height(layers):
        return sum(max(t, FMIN2 if "\n" in n else FMIN) for n, _, t, _, _ in layers)
    tallest = max(height(layers) for _, layers, _ in panels)
    lens = any(layers[0][4] for _, layers, _ in panels)
    below = (mla_depth(FW) if lens else 0.0) + 134.0
    ybase = H / 2.0 + (tallest + TITLE_GAP + 25.0 + FDY - below) / 2.0
    title_y = ybase - tallest - FDY - TITLE_GAP

    prs = Presentation()
    prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H
    SH = prs.slides.add_slide(prs.slide_layouts[6]).shapes
    for x0, layers, title in panels:
        front_device(x0, ybase, layers, title, title_y)
        front_emission(x0, ybase, layers[0][4])
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    prs.save(out)
    print("saved %s  (%d shapes, 0 groups, 0 pictures)" % (out, len(SH)))


figure([(230.0, STACK_A, "Ag reflector + MLA film", 430.0, 340.0),
        (860.0, STACK_B, "DBR + MLA substrate", 1060.0, 330.0)],
       "device_stacks.pptx")

figure([(565.0, STACK_C, "Bottom-emitting OLED", 755.0, 400.0)],
       "device_stack_electrode.pptx")

front_figure([(160.0, STACK_A, "Ag reflector + MLA film"),
              (790.0, STACK_B, "DBR + MLA substrate")],
             "device_stacks_front.pptx")

front_figure([(490.0, STACK_C, "Bottom-emitting OLED")],
             "device_stack_electrode_front.pptx", fw=420.0)
