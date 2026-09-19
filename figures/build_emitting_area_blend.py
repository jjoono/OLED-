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
EMIT_STRENGTH = 26.0                     # design-rule peak emission
WARM = (1.0, 0.92, 0.76)                 # emitted light, like the reference render
SEP = 1.75                               # half distance between the two devices

LAYERS = [  # name, thickness, material spec
    ("Metal cathode",        0.10, dict(base=(0.60, 0.62, 0.65), metallic=1.0, rough=0.32)),
    ("Organic layers",       0.05, dict(base=(0.98, 0.84, 0.60), rough=0.55)),
    ("TCO anode",            0.04, dict(base=(0.74, 0.88, 0.96), rough=0.12, transmission=0.6, ior=1.6)),
    ("Glass substrate",      0.18, dict(base=(0.80, 0.90, 0.98), rough=0.05, transmission=0.9, ior=1.5)),
    ("Outcoupling structure",0.06, dict(base=(0.82, 0.91, 0.98), rough=0.22, transmission=0.9, ior=1.5)),
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

def principled(name, base, metallic=0.0, rough=0.5, transmission=0.0, ior=1.45):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    p = nt.nodes["Principled BSDF"]
    p.inputs["Base Color"].default_value = (*base, 1.0)
    p.inputs["Metallic"].default_value = metallic
    p.inputs["Roughness"].default_value = rough
    p.inputs["IOR"].default_value = ior
    set_in(p, "transmission", transmission)
    if transmission > 0.0:
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

def box(name, sx, sy, sz, cx, cy, z0, mat, c):
    bpy.ops.mesh.primitive_cube_add(size=1.0, location=(cx, cy, z0 + sz / 2.0))
    o = bpy.context.active_object
    o.name = name
    o.scale = (sx, sy, sz)
    bpy.ops.object.transform_apply(scale=True)
    o.data.materials.append(mat)
    link(o, c)
    return o

def look_at(cam, target):
    d = Vector(target) - cam.location
    cam.rotation_euler = d.to_track_quat("-Z", "Y").to_euler()

# ------------------------------------------------------------------ materials
def emitter_material(tag, lam, amp):
    """Emission that falls off as exp(-r / lam) from the panel centre."""
    m = bpy.data.materials.new("Emitting area " + tag)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    tex = nt.nodes.new("ShaderNodeTexCoord")
    ln = nt.nodes.new("ShaderNodeVectorMath"); ln.operation = "LENGTH"
    div = nt.nodes.new("ShaderNodeMath"); div.operation = "DIVIDE"; div.inputs[1].default_value = lam
    neg = nt.nodes.new("ShaderNodeMath"); neg.operation = "MULTIPLY"; neg.inputs[1].default_value = -1.0
    ex = nt.nodes.new("ShaderNodeMath"); ex.operation = "EXPONENT"
    mul = nt.nodes.new("ShaderNodeMath"); mul.operation = "MULTIPLY"; mul.inputs[1].default_value = amp
    em = nt.nodes.new("ShaderNodeEmission")
    em.inputs["Color"].default_value = (*WARM, 1.0)
    tr = nt.nodes.new("ShaderNodeBsdfTransparent")
    mix = nt.nodes.new("ShaderNodeMixShader")
    nt.links.new(tex.outputs["Object"], ln.inputs[0])
    nt.links.new(ln.outputs["Value"], div.inputs[0])
    nt.links.new(div.outputs[0], neg.inputs[0])
    nt.links.new(neg.outputs[0], ex.inputs[0])
    nt.links.new(ex.outputs[0], mul.inputs[0])
    nt.links.new(mul.outputs[0], em.inputs["Strength"])
    nt.links.new(ex.outputs[0], mix.inputs["Fac"])      # transparent where nothing is emitted
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
    dens = N("ShaderNodeMath", operation="MULTIPLY", v1=0.22 * brightness)
    emis = N("ShaderNodeMath", operation="MULTIPLY", v1=2.6 * brightness)
    vol = N("ShaderNodeVolumePrincipled")
    vol.inputs["Color"].default_value = (*WARM, 1.0)
    vol.inputs["Emission Color"].default_value = (*WARM, 1.0)
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
    top = z
    # emitting area: a plane the size of the top face with the exp(-r/lam) falloff
    bpy.ops.mesh.primitive_plane_add(size=1.0, location=(cx, 0.0, top + 0.002))
    e = bpy.context.active_object
    e.name = "Emitting area (%s)" % tag
    e.scale = (PANEL_W, PANEL_D, 1.0)
    bpy.ops.object.transform_apply(scale=True)
    e.data.materials.append(emitter_material(tag, lam, amp))
    e.visible_shadow = False
    link(e, c)
    # light haze above the panel
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
    return top

def studio():
    c = coll("Studio")
    bpy.ops.mesh.primitive_plane_add(size=80.0, location=(0, 0, 0))
    fl = bpy.context.active_object
    fl.name = "Floor"
    fl.data.materials.append(principled("Floor", (0.52, 0.53, 0.55), rough=0.72))
    link(fl, c)
    w = bpy.context.scene.world or bpy.data.worlds.new("World")
    bpy.context.scene.world = w
    w.use_nodes = True
    bg = w.node_tree.nodes["Background"]
    bg.inputs["Color"].default_value = (0.36, 0.37, 0.39, 1.0)
    bg.inputs["Strength"].default_value = 1.0
    for name, loc, energy, size in (("Key light", (-4.0, -6.0, 7.5), 1400.0, 5.0),
                                    ("Fill light", (6.0, -3.0, 5.0), 500.0, 6.0)):
        bpy.ops.object.light_add(type="AREA", location=loc)
        L = bpy.context.active_object
        L.name = name
        L.data.energy = energy
        L.data.size = size
        look_at(L, (0.0, 0.0, 0.3))
        link(L, c)

def cameras():
    c = coll("Cameras")
    cams = {}
    for name, loc, target, lens in (("Cam both",  (0.0, -8.8, 4.4),  (0.0, 0.0, 0.75), 45.0),
                                    ("Cam conventional", (-SEP - 0.2, -5.2, 2.9), (-SEP, 0.0, 0.5), 50.0),
                                    ("Cam design rule",  ( SEP - 0.2, -5.2, 2.9), ( SEP, 0.0, 0.5), 50.0)):
        bpy.ops.object.camera_add(location=loc)
        cam = bpy.context.active_object
        cam.name = name
        cam.data.lens = lens
        look_at(cam, target)
        link(cam, c)
        cams[name] = cam
    return cams

def render_settings(quality):
    s = bpy.context.scene
    s.render.engine = "CYCLES"
    s.cycles.device = "CPU"
    s.cycles.samples = 48 if quality == "test" else 256
    s.cycles.use_adaptive_sampling = True
    s.cycles.adaptive_threshold = 0.02 if quality == "test" else 0.01
    s.cycles.volume_step_rate = 1.0
    s.cycles.volume_max_steps = 256
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
    s.render.resolution_percentage = 50 if quality == "test" else 100
    s.render.image_settings.file_format = "PNG"
    s.render.image_settings.color_mode = "RGB"
    s.view_settings.view_transform = "AgX"
    try:
        s.view_settings.look = "AgX - Medium High Contrast"
    except Exception:
        pass
    s.render.film_transparent = False

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
bpy.context.scene.camera = cams["Cam both"]

blend = os.path.join(OUT, "emitting_area_render.blend")
bpy.ops.wm.save_as_mainfile(filepath=blend, relative_remap=True)
print("saved %s" % blend)
print("spreading %.2fx  ->  lam %.3f vs %.3f   |   light %.2fx  ->  strength %.1f vs %.1f"
      % (SPREAD, LAM_REF / SPREAD, LAM_REF, BRIGHT, EMIT_STRENGTH / BRIGHT, EMIT_STRENGTH))

if DO_RENDER:
    render(cams["Cam both"], 1800, 900, os.path.join(OUT, "emitting_area_render.png"))
    if QUALITY != "test":
        render(cams["Cam conventional"], 1000, 1000,
               os.path.join(OUT, "emitting_area_render_conventional.png"))
        render(cams["Cam design rule"], 1000, 1000,
               os.path.join(OUT, "emitting_area_render_designrule.png"))
