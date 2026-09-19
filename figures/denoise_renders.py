# -*- coding: utf-8 -*-
"""Denoise the Cycles renders with Intel Open Image Denoise (PyPI package `oidn`).

    python3 figures/denoise_renders.py figures/emitting_area_render.png [...]

For every PNG given, the matching multilayer EXR written by build_emitting_area_blend.py
supplies the albedo and normal guide passes; the PNG itself (already through Blender's
AgX view transform) is the colour input, so the denoised result keeps exactly Blender's
look.  Falls back to colour-only denoising when the EXR is missing.  The noisy original
is kept as *_noisy.png.
"""
import os, sys
import numpy as np
from PIL import Image
import oidn

def exr_pass(path, key):
    """Return an (H, W, 3) float32 array for the first layer whose name contains `key`."""
    import OpenEXR, Imath
    f = OpenEXR.InputFile(path)
    hdr = f.header()
    names = [c for c in hdr["channels"] if key in c]
    if len(names) < 3:
        return None
    base = names[0].rsplit(".", 1)[0]
    suff = [".R", ".G", ".B"] if base + ".R" in names else [".X", ".Y", ".Z"]
    dw = hdr["dataWindow"]
    w, h = dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1
    pt = Imath.PixelType(Imath.PixelType.FLOAT)
    chans = [np.frombuffer(f.channel(base + s, pt), dtype=np.float32).reshape(h, w) for s in suff]
    return np.stack(chans, axis=-1).copy()

def _set_flags(flt, **flags):
    """The PyPI wrapper exposes no parameter setter; reach the bundled OIDN 1.x C API."""
    import ctypes, glob
    try:
        so = glob.glob(os.path.join(os.path.dirname(oidn.__file__), "**", "libOpenImageDenoise*.so*"),
                       recursive=True)
        lib = ctypes.CDLL(so[0])
        lib.oidnSetFilter1b.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_bool]
        for k, v in flags.items():
            lib.oidnSetFilter1b(ctypes.c_void_p(flt), k.encode(), bool(v))
    except Exception as e:                      # flags are an optimisation, not a requirement
        print("   (OIDN flags not set: %s)" % e)

def _run(dev, color, w, h, albedo=None, normal=None):
    out = np.zeros_like(color)
    flt = oidn.NewFilter(dev, "RT")
    oidn.SetSharedFilterImage(flt, "color", np.ascontiguousarray(color), oidn.FORMAT_FLOAT3, w, h)
    if albedo is not None:
        oidn.SetSharedFilterImage(flt, "albedo", np.ascontiguousarray(albedo), oidn.FORMAT_FLOAT3, w, h)
    if normal is not None:
        oidn.SetSharedFilterImage(flt, "normal", np.ascontiguousarray(normal), oidn.FORMAT_FLOAT3, w, h)
    oidn.SetSharedFilterImage(flt, "output", out, oidn.FORMAT_FLOAT3, w, h)
    _set_flags(flt, hdr=False, srgb=True)
    oidn.CommitFilter(flt)
    oidn.ExecuteFilter(flt)
    oidn.ReleaseFilter(flt)
    return out

def denoise(png):
    src = Image.open(png)
    alpha = np.asarray(src.split()[-1], dtype=np.float32) / 255.0 if src.mode in ("RGBA", "LA") else None
    img = np.asarray(src.convert("RGB"), dtype=np.float32) / 255.0
    h, w, _ = img.shape
    exr = os.path.splitext(png)[0] + ".exr"
    albedo = normal = None
    if os.path.exists(exr):
        albedo = exr_pass(exr, "Denoising Albedo")
        normal = exr_pass(exr, "Denoising Normal")
        if normal is not None:
            normal = np.clip(normal, -1.0, 1.0)
    dev = oidn.NewDevice()
    oidn.CommitDevice(dev)
    out = _run(dev, img, w, h, albedo, normal)
    if alpha is not None:
        # the glow's noise lives in alpha once the compositor rebuilds it there,
        # so that channel needs denoising too (as a grey triplet, no guides)
        a3 = _run(dev, np.repeat(alpha[:, :, None], 3, axis=2), w, h)
        alpha = np.clip(a3[:, :, 0], 0.0, 1.0)
    err = oidn.GetDeviceError(dev)
    oidn.ReleaseDevice(dev)
    if err and err[0] != oidn.ERROR_NONE:
        raise RuntimeError("OIDN: %s" % (err,))
    noisy = os.path.splitext(png)[0] + "_noisy.png"
    os.replace(png, noisy)
    res = Image.fromarray((np.clip(out, 0, 1) * 255 + 0.5).astype(np.uint8))
    if alpha is not None:
        res.putalpha(Image.fromarray((alpha * 255 + 0.5).astype(np.uint8)))
    res.save(png)
    print("denoised %s   guides: albedo=%s normal=%s   alpha=%s" % (
        os.path.basename(png), albedo is not None, normal is not None, alpha is not None))

if __name__ == "__main__":
    for p in sys.argv[1:]:
        denoise(p)
