# -*- coding: utf-8 -*-
"""Build (and optionally render) the emitting-area comparison as a Blender scene.

    blender -b -P figures/build_emitting_area_blend.py -- --out figures [--quality test|final]

Two OLED stacks on a grey studio floor.  The lossy device (left) glows in a
small dim spot under a short haze cone; the design-rule device (right) glows
across the whole panel under a tall bright cone.  Nothing about the two glows
is eyeballed: the same round-trip loss model as Fig.1(a) supplies

    spreading length ratio  2.32x  ->  radius of the emission falloff
    total light ratio       1.84x  ->  emission strength

Cycles on CPU, OpenImageDenoise if the build has it.  Saves emitting_area_render.blend
plus three renders (both devices, and one close-up per device).
"""
import bpy, math, sys, os, time
from mathutils import Vector

# ------------------------------------------------------------------ arguments
argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
OUT = os.path.abspath(argv[argv.index("--out") + 1]) if "--out" in argv else os.getcwd()
QUALITY = argv[argv.index("--quality") + 1] if "--quality" in argv else "final"
DO_RENDER = "--no-render" not in argv
BG = argv[argv.index("--bg") + 1] if "--bg" in argv else "transparent"
ONLY_CAM = argv[argv.index("--cam") + 1] if "--cam" in argv else None
RES_PCT = int(argv[argv.index("--res") + 1]) if "--res" in argv else None
SAMPLES = int(argv[argv.index("--samples") + 1]) if "--samples" in argv else None
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ loss model
def model(eta, r_met, t_tco):
    rho = (1 - eta) * r_met * t_tco ** 2
    return 1.0 / math.log(1.0 / rho), eta / (1.0 - rho)

LAM_C, TOT_C = model(.30, .72, .90)      # conventional
LAM_D, TOT_D = model(.30, .98, .995)     # design rule
SPREAD, BRIGHT = LAM_D / LAM_C, TOT_D / TOT_C

PANEL_W, PANEL_D = 2.0, 1.5              # device footprint, Blender units
LAM_REF = 0.60                           # design-rule spreading length (0.30 panel widths)
EMIT_STRENGTH = 9.0
LENS_EMIT_FRAC = 0.85     # how much of the glow the caps themselves carry                     # design-rule peak emission
GLOW = (0.06, 1.0, 0.30)                 # emitted light: vivid green (linear)
HAZE_SCATTER = (0.55, 1.0, 0.72)         # what the haze scatters, kept lighter
SEP = 1.75                               # half distance between the two devices
LENS_PITCH = 0.125                       # micro-lens centre-to-centre, hexagonal packing
LENS_FILL = 0.97                         # lens radius as a fraction of half the pitch
LENS_SINK = 0.16                         # how far each cap sits in the film, x radius
LENS_SEGS, LENS_RINGS = 28, 9

LAYERS = [  # name, thickness, material spec
    ("Metal cathode",        0.10, dict(base=(0.60, 0.62, 0.65), metallic=1.0, rough=0.32)),
    ("Organic layers",       0.05, dict(base=(0.98, 0.84, 0.60), rough=0.55)),
    ("TCO anode",            0.04, dict(base=(0.74, 0.88, 0.96), rough=0.12, transmission=0.6, ior=1.6)),
    ("Glass substrate",      0.18, dict(base=(0.80, 0.90, 0.98), rough=0.05, transmission=0.6, ior=1.5)),
    ("Outcoupling film",     0.07, dict(base=(0.94, 0.97, 1.00), rough=0.12, transmission=0.10,
                                        ior=1.5)),
]

# ------------------------------------------------------------------ helpers
def clear():
    bpy.ops.wm.read_factory_settings(use_empty=True)

def coll(name, parent=None):
    c = bpy.data.collections.new(name)
    (parent or bpy.context.scene.collection).children.link(c)
    return c

def link(obj, c):
    for old in obj.users_collection:
        old.objects.unlink(obj)
    c.objects.link(obj)

def set_in(node, name, value):
    """Set a node input by name, tolerating 4.0 renames."""
    aliases = {"transmission": ("Transmission Weight", "Transmission"),
               "emission_color": ("Emission Color", "Emission"),
               "emission_strength": ("Emission Strength",),
               "specular": ("Specular IOR Level", "Specular")}
    for key in aliases.get(name, (name,)):
        if key in node.inputs:
            node.inputs[key].default_value = value
            return
    raise KeyError(name)

def dome(r, segs, rings):
    """Vertices and faces of one spherical cap, pole up, closed underneath."""
    verts = [(0.0, 0.0, r)]
    for i in range(1, rings + 1):
        th = (math.pi / 2.0) * i / rings
        st, ct = math.sin(th), math.cos(th)
        for j in range(segs):
            ph = 2.0 * math.pi * j / segs
            verts.append((r * st * math.cos(ph), r * st * math.sin(ph), r * ct))
    faces = [(0, 1 + j, 1 + (j + 1) % segs) for j in range(segs)]
    for i in range(1, rings):
        a0, b0 = 1 + (i - 1) * segs, 1 + i * segs
        faces += [(a0 + j, b0 + j, b0 + (j + 1) % segs, a0 + (j + 1) % segs) for j in range(segs)]
    rim = [1 + (rings - 1) * segs + j for j in range(segs)]
    faces.append(tuple(reversed(rim)))
    return verts, faces


def lens_array(cx, z_top, tag, c, mat):
    """Micro-lens array: spherical caps on a hexagonal lattice, one joined mesh.

    Rows are offset by half a pitch and spaced by pitch*sqrt(3)/2, which is the
    hexagonal (triangular) packing an MLA actually uses.  Each cap is sunk a little
    into the film so the rim is hidden and the lens reads as a cap, not a ball.
    """
    r = LENS_PITCH * 0.5 * LENS_FILL
    row_h = LENS_PITCH * math.sqrt(3.0) / 2.0
    dv, df = dome(r, LENS_SEGS, LENS_RINGS)
    z0 = z_top - LENS_SINK * r
    verts, faces, n = [], [], 0
    jmax = int(((PANEL_D / 2.0) - r) / row_h) + 1
    imax = int(((PANEL_W / 2.0) - r) / LENS_PITCH) + 1
    for j in range(-jmax, jmax + 1):
        y = j * row_h
        if abs(y) > PANEL_D / 2.0 - r:
            continue
        xoff = (LENS_PITCH / 2.0) if j % 2 else 0.0
        for i in range(-imax, imax + 1):
            x = i * LENS_PITCH + xoff
            if abs(x) > PANEL_W / 2.0 - r:
                continue
            verts += [(vx + x, vy + y, vz + z0) for vx, vy, vz in dv]
            faces += [tuple(k + n for k in f) for f in df]
            n += len(dv)
    me = bpy.data.meshes.new("Micro-lens array (%s)" % tag)
    me.from_pydata(verts, [], faces)
    me.update()
    for poly in me.polygons:
        poly.use_smooth = True
    me.materials.append(mat)
    o = bpy.data.objects.new("Micro-lens array (%s)" % tag, me)
    bpy.context.scene.collection.objects.link(o)
    o.location.x = cx
    link(o, c)
    return o, z_top + r * (1.0 - LENS_SINK), n // len(dv)


def principled(name, base, metallic=0.0, rough=0.5, transmission=0.0, ior=1.45,
               shadow_through=True):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    p = nt.nodes["Principled BSDF"]
    p.inputs["Base Color"].default_value = (*base, 1.0)
    p.inputs["Metallic"].default_value = metallic
    p.inputs["Roughness"].default_value = rough
    p.inputs["IOR"].default_value = ior
    set_in(p, "transmission", transmission)
    if transmission > 0.0 and shadow_through:
        # shadow rays go straight through, so light reaches the layers underneath
        out = nt.nodes["Material Output"]
        lp = nt.nodes.new("ShaderNodeLightPath")
        tr = nt.nodes.new("ShaderNodeBsdfTransparent")
        mix = nt.nodes.new("ShaderNodeMixShader")
        nt.links.new(lp.outputs["Is Shadow Ray"], mix.inputs["Fac"])
        nt.links.new(p.outputs[0], mix.inputs[1])
        nt.links.new(tr.outputs[0], mix.inputs[2])
        nt.links.new(mix.outputs[0], out.inputs["Surface"])
    return m

def box(name, sx, sy, sz, cx, cy, z0, mat, c, bevel=0.0035):
    bpy.ops.mesh.primitive_cube_add(size=1.0, location=(cx, cy, z0 + sz / 2.0))
    o = bpy.context.active_object
    o.name = name
    o.scale = (sx, sy, sz)
    bpy.ops.object.transform_apply(scale=True)
    o.data.materials.append(mat)
    if bevel:
        m = o.modifiers.new("Bevel", "BEVEL")
        m.width = min(bevel, sz * 0.22)      # thin layers get a proportionally smaller one
        m.segments = 4
        m.limit_method = "ANGLE"
        m.angle_limit = math.radians(30.0)
        m.use_clamp_overlap = True
    link(o, c)
    return o

def look_at(cam, target):
    d = Vector(target) - cam.location
    cam.rotation_euler = d.to_track_quat("-Z", "Y").to_euler()

# ------------------------------------------------------------------ materials
def radial_falloff(nt, lam, amp):
    """exp(-r/lam) * amp from the panel centre, in the object's own XY plane."""
    tex = nt.nodes.new("ShaderNodeTexCoord")
    sep = nt.nodes.new("ShaderNodeSeparateXYZ")
    flat = nt.nodes.new("ShaderNodeCombineXYZ")
    ln = nt.nodes.new("ShaderNodeVectorMath"); ln.operation = "LENGTH"
    div = nt.nodes.new("ShaderNodeMath"); div.operation = "DIVIDE"; div.inputs[1].default_value = lam
    neg = nt.nodes.new("ShaderNodeMath"); neg.operation = "MULTIPLY"; neg.inputs[1].default_value = -1.0
    ex = nt.nodes.new("ShaderNodeMath"); ex.operation = "EXPONENT"
    mul = nt.nodes.new("ShaderNodeMath"); mul.operation = "MULTIPLY"; mul.inputs[1].default_value = amp
    L = nt.links.new
    L(tex.outputs["Object"], sep.inputs[0])
    L(sep.outputs["X"], flat.inputs["X"]); L(sep.outputs["Y"], flat.inputs["Y"])
    L(flat.outputs[0], ln.inputs[0])
    L(ln.outputs["Value"], div.inputs[0])
    L(div.outputs[0], neg.inputs[0]); L(neg.outputs[0], ex.inputs[0]); L(ex.outputs[0], mul.inputs[0])
    return ex.outputs[0], mul.outputs[0]        # (0..1 shape, shape * amp)


def lens_material(tag, lam, amp, strength):
    """White cap that also emits: light really does leave through the lenses, so
    they can stay opaque enough to shade properly and still carry the glow."""
    m = principled("Micro-lens %s" % tag, base=(0.95, 0.97, 1.00), rough=0.12,
                   transmission=0.08, ior=1.50, shadow_through=False)
    nt = m.node_tree
    p = nt.nodes["Principled BSDF"]
    shape, scaled = radial_falloff(nt, lam, amp * strength)
    set_in(p, "emission_color", (*GLOW, 1.0))
    nt.links.new(scaled, p.inputs["Emission Strength"])
    return m


def emitter_material(tag, lam, amp):
    """Emission that falls off as exp(-r / lam) from the panel centre."""
    m = bpy.data.materials.new("Emitting area " + tag)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    shape, scaled = radial_falloff(nt, lam, amp)
    em = nt.nodes.new("ShaderNodeEmission")
    em.inputs["Color"].default_value = (*GLOW, 1.0)
    tr = nt.nodes.new("ShaderNodeBsdfTransparent")
    mix = nt.nodes.new("ShaderNodeMixShader")
    nt.links.new(scaled, em.inputs["Strength"])
    nt.links.new(shape, mix.inputs["Fac"])             # transparent where nothing is emitted
    nt.links.new(tr.outputs[0], mix.inputs[1])
    nt.links.new(em.outputs[0], mix.inputs[2])
    nt.links.new(mix.outputs[0], out.inputs["Surface"])
    m.blend_method = "BLEND"
    return m

def haze_material(tag, height, r_base, slope, brightness):
    """Thin glowing volume: fades with height and softens toward the cone wall."""
    m = bpy.data.materials.new("Light haze " + tag)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    N = lambda kind, **kw: _node(nt, kind, **kw)
    out = N("ShaderNodeOutputMaterial")
    tex = N("ShaderNodeTexCoord")
    sep = N("ShaderNodeSeparateXYZ")
    # vertical fade  (1 - z/h)^2.2
    zn = N("ShaderNodeMath", operation="DIVIDE", v1=height)
    inv = N("ShaderNodeMath", operation="SUBTRACT", v0=1.0)
    vfade = N("ShaderNodeMath", operation="POWER", v1=2.2)
    # radial fade  1 - (r / R(z))^3   with  R(z) = r_base + slope * z
    xy = N("ShaderNodeCombineXYZ")
    rad = N("ShaderNodeVectorMath", operation="LENGTH")
    rz = N("ShaderNodeMath", operation="MULTIPLY_ADD", v1=slope, v2=r_base)
    t = N("ShaderNodeMath", operation="DIVIDE")
    t3 = N("ShaderNodeMath", operation="POWER", v1=3.0)
    rfade = N("ShaderNodeMath", operation="SUBTRACT", v0=1.0, clamp=True)
    fade = N("ShaderNodeMath", operation="MULTIPLY")
    dens = N("ShaderNodeMath", operation="MULTIPLY", v1=0.26 * brightness)
    emis = N("ShaderNodeMath", operation="MULTIPLY", v1=0.85 * brightness)
    vol = N("ShaderNodeVolumePrincipled")
    vol.inputs["Color"].default_value = (*HAZE_SCATTER, 1.0)
    vol.inputs["Emission Color"].default_value = (*GLOW, 1.0)
    vol.inputs["Anisotropy"].default_value = 0.3
    L = nt.links.new
    L(tex.outputs["Object"], sep.inputs[0])
    L(sep.outputs["Z"], zn.inputs[0]); L(zn.outputs[0], inv.inputs[1]); L(inv.outputs[0], vfade.inputs[0])
    L(sep.outputs["X"], xy.inputs["X"]); L(sep.outputs["Y"], xy.inputs["Y"])
    L(xy.outputs[0], rad.inputs[0])
    L(sep.outputs["Z"], rz.inputs[0])
    L(rad.outputs["Value"], t.inputs[0]); L(rz.outputs[0], t.inputs[1])
    L(t.outputs[0], t3.inputs[0]); L(t3.outputs[0], rfade.inputs[1])
    L(vfade.outputs[0], fade.inputs[0]); L(rfade.outputs[0], fade.inputs[1])
    L(fade.outputs[0], dens.inputs[0]); L(fade.outputs[0], emis.inputs[0])
    L(dens.outputs[0], vol.inputs["Density"]); L(emis.outputs[0], vol.inputs["Emission Strength"])
    L(vol.outputs[0], out.inputs["Volume"])
    return m

def _node(nt, kind, operation=None, v0=None, v1=None, v2=None, clamp=False):
    n = nt.nodes.new(kind)
    if operation:
        n.operation = operation
    for i, v in enumerate((v0, v1, v2)):
        if v is not None:
            n.inputs[i].default_value = v
    if clamp:
        n.use_clamp = True
    return n

# ------------------------------------------------------------------ scene
def device(cx, tag, lam, amp, cone_r, cone_h, brightness):
    c = coll("Device " + tag)
    z = 0.0
    for name, t, spec in LAYERS:
        box("%s (%s)" % (name, tag), PANEL_W, PANEL_D, t, cx, 0.0, z,
            principled("%s %s" % (name, tag), **spec), c)
        z += t
    film = z                                   # top of the outcoupling film
    lens_mat = lens_material(tag, lam, amp, LENS_EMIT_FRAC)
    _, top, n_lens = lens_array(cx, film, tag, c, lens_mat)
    # emitting area: sits under the lens array, so the light leaves through the lenses
    bpy.ops.mesh.primitive_plane_add(size=1.0, location=(cx, 0.0, film + 0.0015))
    e = bpy.context.active_object
    e.name = "Emitting area (%s)" % tag
    e.scale = (PANEL_W, PANEL_D, 1.0)
    bpy.ops.object.transform_apply(scale=True)
    e.data.materials.append(emitter_material(tag, lam, amp))
    e.visible_shadow = False
    link(e, c)
    # light haze above the lenses
    slope = 0.55
    bpy.ops.mesh.primitive_cone_add(vertices=72, radius1=cone_r, radius2=cone_r + cone_h * slope,
                                    depth=cone_h, end_fill_type="NGON",
                                    location=(cx, 0.0, top + cone_h / 2.0))
    k = bpy.context.active_object
    k.name = "Light haze (%s)" % tag
    for v in k.data.vertices:                     # origin at the base so Object Z runs 0..h
        v.co.z += cone_h / 2.0
    k.location.z = top
    k.data.materials.append(haze_material(tag, cone_h, cone_r, slope, brightness))
    k.visible_shadow = False
    k.visible_diffuse = False
    k.visible_glossy = False
    link(k, c)
    print("   %-12s  %d micro-lenses, pitch %.3f" % (tag, n_lens, LENS_PITCH))
    return top


def studio():
    c = coll("Studio")
    bpy.ops.mesh.primitive_plane_add(size=80.0, location=(0, 0, 0))
    fl = bpy.context.active_object
    fl.name = "Floor"
    fl.data.materials.append(principled("Floor", (0.52, 0.53, 0.55), rough=0.72))
    if BG == "transparent":
        fl.is_shadow_catcher = True      # keeps the contact shadow, drops the grey plane
    link(fl, c)

    # a graded environment rather than a flat grey: the lens caps and the glass
    # pick up a rolloff instead of one constant value, which is most of the
    # difference between a render that looks shot and one that looks computed
    w = bpy.context.scene.world or bpy.data.worlds.new("World")
    bpy.context.scene.world = w
    w.use_nodes = True
    nt = w.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputWorld")
    bg = nt.nodes.new("ShaderNodeBackground")
    bg.inputs["Strength"].default_value = 0.11
    tex = nt.nodes.new("ShaderNodeTexCoord")
    sep = nt.nodes.new("ShaderNodeSeparateXYZ")
    mr = nt.nodes.new("ShaderNodeMapRange")
    mr.inputs["From Min"].default_value = -0.9
    mr.inputs["From Max"].default_value = 0.9
    ramp = nt.nodes.new("ShaderNodeValToRGB")
    ramp.color_ramp.interpolation = "B_SPLINE"
    for pos, col in ((0.0, (0.16, 0.17, 0.19)), (0.45, (0.42, 0.44, 0.48)),
                     (1.0, (0.92, 0.94, 0.97))):
        e = ramp.color_ramp.elements.new(pos) if pos not in (0.0,) else ramp.color_ramp.elements[0]
        e.position, e.color = pos, (*col, 1.0)
    nt.links.new(tex.outputs["Generated"], sep.inputs[0])
    nt.links.new(sep.outputs["X"], mr.inputs["Value"])   # sideways, not overhead
    nt.links.new(mr.outputs[0], ramp.inputs["Fac"])
    nt.links.new(ramp.outputs["Color"], bg.inputs["Color"])
    nt.links.new(bg.outputs[0], out.inputs["Surface"])

    # big soft key, wide low fill, and a rim behind to separate the slab edges
    for name, loc, energy, size, target in (
            ("Key light",  (-2.6, -3.6, 12.0), 1500.0, 9.0, (0.0, 0.0, 0.45)),
            ("Fill light", ( 6.5, -5.5, 3.4),  420.0, 10.0, (0.0, 0.0, 0.32)),
            ("Rim light",  ( 1.5,  7.5, 3.0),  260.0, 6.0, (0.0, 0.0, 0.40))):
        bpy.ops.object.light_add(type="AREA", location=loc)
        L = bpy.context.active_object
        L.name = name
        L.data.energy = energy
        L.data.size = size
        L.data.spread = math.radians(120.0)
        if name != "Key light":               # a single shadow keeps the floor clean
            for owner, attr in ((L.data, "use_shadow"), (getattr(L.data, "cycles", None), "cast_shadow")):
                if owner is not None and hasattr(owner, attr):
                    setattr(owner, attr, False)
                    break
        look_at(L, target)
        link(L, c)


def cameras():
    c = coll("Cameras")
    cams = {}
    for name, loc, target, lens in (("Cam both",  (0.0, -15.40, 9.20), (0.0, 0.0, 0.70), 85.0),
                                    ("Cam conventional", (-SEP - 0.34, -7.83, 5.66), (-SEP, 0.0, 0.42), 85.0),
                                    ("Cam design rule",  ( SEP - 0.34, -7.83, 5.66), ( SEP, 0.0, 0.42), 85.0)):
        bpy.ops.object.camera_add(location=loc)
        cam = bpy.context.active_object
        cam.name = name
        cam.data.lens = lens
        look_at(cam, target)
        link(cam, c)
        cams[name] = cam
    return cams

def compositor():
    """Give the emissive haze an alpha so the glow survives a transparent background.

    Cycles gives a thin emissive volume almost no alpha, so with film_transparent the
    cone is there in RGB but vanishes once composited.  Alpha is rebuilt from the
    per-pixel channel maximum -- not the luminance, which for a saturated colour would
    push the strong channel past 1 on un-premultiply and clip the hue back to white.
    The colour is un-premultiplied only where the render had no alpha of its own, so
    the glow composites cleanly over a light page and over a dark one alike.
    """
    s = bpy.context.scene
    s.use_nodes = True
    nt = s.node_tree
    nt.nodes.clear()
    rl = nt.nodes.new("CompositorNodeRLayers")
    sat = nt.nodes.new("CompositorNodeHueSat")       # OCIO "looks" are unavailable in
    for key, val in (("Saturation", 1.30), ("Fac", 1.0)):   # this build, so do it here
        if key in sat.inputs:
            sat.inputs[key].default_value = val
    try:
        sep = nt.nodes.new("CompositorNodeSeparateColor")
        sep.mode = "RGB"
        chans = ("Red", "Green", "Blue")
    except RuntimeError:
        sep = nt.nodes.new("CompositorNodeSepRGBA")
        chans = ("R", "G", "B")
    m1 = nt.nodes.new("CompositorNodeMath"); m1.operation = "MAXIMUM"
    lum = nt.nodes.new("CompositorNodeMath"); lum.operation = "MAXIMUM"
    safe = nt.nodes.new("CompositorNodeMath"); safe.operation = "MAXIMUM"
    safe.inputs[1].default_value = 0.004
    unp = nt.nodes.new("CompositorNodeMixRGB"); unp.blend_type = "DIVIDE"
    unp.inputs[0].default_value = 1.0
    hard = nt.nodes.new("CompositorNodeMath")        # partial alpha (shadow) -> 1
    hard.operation = "MULTIPLY"; hard.inputs[1].default_value = 12.0; hard.use_clamp = True
    keep = nt.nodes.new("CompositorNodeMixRGB"); keep.blend_type = "MIX"
    amax = nt.nodes.new("CompositorNodeMath"); amax.operation = "MAXIMUM"
    sa = nt.nodes.new("CompositorNodeSetAlpha"); sa.mode = "REPLACE_ALPHA"
    out = nt.nodes.new("CompositorNodeComposite")
    L = nt.links.new
    L(rl.outputs["Image"], sat.inputs["Image"])
    L(sat.outputs[0], sep.inputs[0])
    L(sep.outputs[chans[0]], m1.inputs[0])
    L(sep.outputs[chans[1]], m1.inputs[1])
    L(m1.outputs[0], lum.inputs[0])
    L(sep.outputs[chans[2]], lum.inputs[1])
    L(lum.outputs[0], safe.inputs[0])
    L(sat.outputs[0], unp.inputs[1])
    L(safe.outputs[0], unp.inputs[2])
    L(rl.outputs["Alpha"], hard.inputs[0])
    L(hard.outputs[0], keep.inputs[0])
    L(unp.outputs[0], keep.inputs[1])
    L(sat.outputs[0], keep.inputs[2])
    L(lum.outputs[0], amax.inputs[0])
    L(rl.outputs["Alpha"], amax.inputs[1])
    L(keep.outputs[0], sa.inputs["Image"])
    L(amax.outputs[0], sa.inputs["Alpha"])
    L(sa.outputs["Image"], out.inputs["Image"])
    for i, n in enumerate(nt.nodes):
        n.location = (260 * i, 0)

def render_settings(quality):
    s = bpy.context.scene
    s.render.engine = "CYCLES"
    s.cycles.device = "CPU"
    s.cycles.samples = SAMPLES or (64 if quality == "test" else 320)
    s.cycles.use_adaptive_sampling = True
    s.cycles.adaptive_threshold = 0.02 if quality == "test" else 0.006
    s.cycles.volume_step_rate = 1.0 if quality == "test" else 0.5
    s.cycles.volume_max_steps = 256 if quality == "test" else 1024
    # the lens caps sit on glass on glass: give transmission room, but keep
    # caustics off and clamp indirect so the image stays free of fireflies
    s.cycles.max_bounces = 16
    s.cycles.transmission_bounces = 16
    s.cycles.glossy_bounces = 8
    s.cycles.diffuse_bounces = 4
    s.cycles.transparent_max_bounces = 24
    s.cycles.caustics_reflective = False
    s.cycles.caustics_refractive = False
    s.cycles.blur_glossy = 0.0
    s.cycles.sample_clamp_indirect = 24.0
    for owner in (s.cycles, s.render):    # the pixel filter moved between versions
        if hasattr(owner, "filter_width"):
            owner.filter_width = 1.20     # a touch crisper than the 1.5 default
            break
    s.render.use_persistent_data = True
    have = [i.identifier for i in s.cycles.bl_rna.properties["denoiser"].enum_items]
    if "OPENIMAGEDENOISE" in have:
        s.cycles.use_denoising = True
        s.cycles.denoiser = "OPENIMAGEDENOISE"
        s.cycles.denoising_use_gpu = False
    else:                                   # e.g. the Ubuntu apt build: no OIDN
        s.cycles.use_denoising = False      # -> denoise_renders.py does it afterwards
    bpy.context.view_layer.cycles.denoising_store_passes = True   # albedo + normal for OIDN
    print("   denoiser: %s   samples: %d" % ("OIDN" if s.cycles.use_denoising else "none",
                                               s.cycles.samples))
    s.render.resolution_percentage = RES_PCT or (40 if quality == "test" else 100)
    s.render.image_settings.file_format = "PNG"
    s.render.image_settings.color_mode = "RGBA" if BG == "transparent" else "RGB"
    s.render.image_settings.compression = 20
    s.view_settings.view_transform = "AgX"
    try:
        s.view_settings.look = "AgX - Medium High Contrast"
    except Exception:
        pass
    s.render.film_transparent = (BG == "transparent")


def render(cam, w, h, path):
    s = bpy.context.scene
    s.camera = cam
    s.render.resolution_x, s.render.resolution_y = w, h
    s.render.filepath = path
    t = time.time()
    bpy.ops.render.render(write_still=True)
    # same result again as a linear multilayer EXR carrying the denoising passes
    fmt, depth = s.render.image_settings.file_format, s.render.image_settings.color_depth
    s.render.image_settings.file_format = "OPEN_EXR_MULTILAYER"
    s.render.image_settings.color_depth = "32"
    bpy.data.images["Render Result"].save_render(os.path.splitext(path)[0] + ".exr", scene=s)
    s.render.image_settings.file_format, s.render.image_settings.color_depth = fmt, depth
    print("   rendered %s  (%.0f s)" % (os.path.basename(path), time.time() - t))

# ------------------------------------------------------------------ build
clear()
studio()
device(-SEP, "conventional", LAM_REF / SPREAD, EMIT_STRENGTH / BRIGHT,
       cone_r=0.48, cone_h=1.25, brightness=1.0 / BRIGHT)
device(+SEP, "design rule", LAM_REF, EMIT_STRENGTH,
       cone_r=0.92, cone_h=2.05, brightness=1.0)
cams = cameras()
render_settings(QUALITY)
if BG == "transparent":
    compositor()
bpy.context.scene.camera = cams["Cam both"]

blend = os.path.join(OUT, "emitting_area_render.blend")
bpy.ops.wm.save_as_mainfile(filepath=blend, relative_remap=True)
print("saved %s" % blend)
print("spreading %.2fx  ->  lam %.3f vs %.3f   |   light %.2fx  ->  strength %.1f vs %.1f"
      % (SPREAD, LAM_REF / SPREAD, LAM_REF, BRIGHT, EMIT_STRENGTH / BRIGHT, EMIT_STRENGTH))

# 2400 px wide is a Nature double-column figure (180 mm) at 300 dpi
JOBS = [("Cam both", 2400, 1200, "emitting_area_render.png"),
        ("Cam conventional", 1200, 1200, "emitting_area_render_conventional.png"),
        ("Cam design rule", 1200, 1200, "emitting_area_render_designrule.png")]
if DO_RENDER:
    for cam_name, w, h, fn in JOBS:
        if ONLY_CAM and ONLY_CAM.lower() not in cam_name.lower():
            continue
        if not ONLY_CAM and QUALITY == "test" and cam_name != "Cam both":
            continue
        render(cams[cam_name], w, h, os.path.join(OUT, fn))
