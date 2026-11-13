"""
apron2.py -- "clear-core / scattering-apron" variant of apron.py: the SL under the electrically active
square is (nearly) non-scattering, so obliquely emitted light walks laterally out of the lossy active
zone in a few bounces; the SL under the electrode-free apron is scattering and randomises it there.
Zone-dependent scattering coefficient, SL thickness d.  Exponential free paths are re-sampled when a
photon crosses the active-zone boundary (memoryless, so this is exact).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recycling as rc  # noqa: E402
from devices import run_proposed  # noqa: E402
import optimize as opt  # noqa: E402
import analyze as an  # noqa: E402
from apron import ce_from_arrivals, apron_tables  # noqa: E402


def _box_dist(x, y, ux, uy, h):
    """distance along (ux,uy) to the boundary of the square |x|,|y|<=h (exit if inside, entry if outside)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t1x = (-h - x) / ux; t2x = (h - x) / ux
        t1y = (-h - y) / uy; t2y = (h - y) / uy
    tminx = np.where(ux != 0, np.minimum(t1x, t2x), np.where(np.abs(x) <= h, -np.inf, np.inf))
    tmaxx = np.where(ux != 0, np.maximum(t1x, t2x), np.where(np.abs(x) <= h, np.inf, -np.inf))
    tminy = np.where(uy != 0, np.minimum(t1y, t2y), np.where(np.abs(y) <= h, -np.inf, np.inf))
    tmaxy = np.where(uy != 0, np.maximum(t1y, t2y), np.where(np.abs(y) <= h, np.inf, -np.inf))
    tn = np.maximum(tminx, tminy); tf = np.minimum(tmaxx, tmaxy)
    inside = (np.abs(x) <= h) & (np.abs(y) <= h)
    out = np.where(inside, tf, np.where((tn <= tf) & (tn > 1e-9), tn, np.inf))
    return np.where(out > 1e-9, out, np.inf)


def mc_cell2(tA, tP, th_deg, wl, Pin, qw, d_um, musA, musP, g, L_um, a_um, R_wall=0.95,
             n_photons=60000, max_steps=60000, seed=5):
    rng = np.random.default_rng(seed)
    nwl, nth = len(wl), len(th_deg)
    W = Pin * qw[None, :, None]
    total_in = W.sum()
    cdf = np.cumsum(W.ravel()); cdf /= cdf[-1]
    pol, il, it = np.unravel_index(np.searchsorted(cdf, rng.random(n_photons)), W.shape)
    dth = th_deg[1] - th_deg[0]
    theta = np.deg2rad(th_deg[it] + (rng.random(n_photons) - 0.5) * dth)
    phi = 2 * np.pi * rng.random(n_photons)
    ux, uy, uz = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), -np.cos(theta)
    h = a_um / 2; half = L_um / 2
    x = (rng.random(n_photons) - 0.5) * a_um * (1 - 1e-9)
    y = (rng.random(n_photons) - 0.5) * a_um * (1 - 1e-9)
    z = np.full(n_photons, d_um)
    alive = np.ones(n_photons, bool)
    edges = np.deg2rad(np.linspace(0, 90, nth + 1))
    arrA = np.zeros((2, nwl, nth)); arrP = np.zeros((2, nwl, nth))
    fate = dict(escA=0, escP=0, topabsA=0, topabsP=0, bot=0, wall=0)
    for step in range(max_steps):
        ids = np.nonzero(alive)[0]
        if ids.size == 0:
            break
        xi_, yi_, uxi, uyi, uzi = x[ids], y[ids], ux[ids], uy[ids], uz[ids]
        act = (np.abs(xi_) <= h) & (np.abs(yi_) <= h)
        mus = np.where(act, musA, musP)
        dz = np.where(uzi > 0, (d_um - z[ids]) / np.maximum(uzi, 1e-300), z[ids] / np.maximum(-uzi, 1e-300))
        with np.errstate(divide="ignore", invalid="ignore"):
            dwx = np.where(uxi > 0, (half - xi_) / uxi, np.where(uxi < 0, (-half - xi_) / uxi, np.inf))
            dwy = np.where(uyi > 0, (half - yi_) / uyi, np.where(uyi < 0, (-half - yi_) / uyi, np.inf))
        dbx = _box_dist(xi_, yi_, uxi, uyi, h)
        s = -np.log(rng.random(ids.size)) / mus
        dmin = np.minimum.reduce([dz, dwx, dwy, dbx, s])
        x[ids] += uxi * dmin; y[ids] += uyi * dmin; z[ids] += uzi * dmin
        kind = np.select([dmin == dz, (dmin == dwx) | (dmin == dwy), dmin == dbx], [0, 1, 2], 3)
        # zone boundary: nudge across and continue (free path re-sampled next step)
        bb = ids[kind == 2]
        x[bb] += ux[bb] * 1e-6; y[bb] += uy[bb] * 1e-6
        # walls
        wids = ids[kind == 1]
        if wids.size:
            ok = rng.random(wids.size) < R_wall
            fate["wall"] += int((~ok).sum()); alive[wids[~ok]] = False
            w2 = wids[ok]
            hx = np.abs(np.abs(x[w2]) - half) < 1e-7
            ux[w2[hx]] *= -1; uy[w2[~hx]] *= -1
            x[w2] = np.clip(x[w2], -half + 1e-7, half - 1e-7); y[w2] = np.clip(y[w2], -half + 1e-7, half - 1e-7)
        bz = ids[kind == 0]
        if bz.size:
            up = uz[bz] > 0
            ht = bz[up]
            if ht.size:
                z[ht] = d_um
                tb = np.minimum(np.searchsorted(edges, np.arccos(np.clip(uz[ht], 0, 1))) - 1, nth - 1)
                a_ = (np.abs(x[ht]) <= h) & (np.abs(y[ht]) <= h)
                Rv = np.where(a_, tA["Rt"][pol[ht], il[ht], tb], tP["Rt"][pol[ht], il[ht], tb])
                Tv = np.where(a_, tA["Tt"][pol[ht], il[ht], tb], tP["Tt"][pol[ht], il[ht], tb])
                np.add.at(arrA, (pol[ht[a_]], il[ht[a_]], tb[a_]), 1.0)
                np.add.at(arrP, (pol[ht[~a_]], il[ht[~a_]], tb[~a_]), 1.0)
                r = rng.random(ht.size)
                refl = r < Rv; trans = (~refl) & (r < Rv + Tv)
                uz[ht[refl]] *= -1
                fate["escA"] += int((trans & a_).sum()); fate["escP"] += int((trans & ~a_).sum())
                fate["topabsA"] += int((~refl & ~trans & a_).sum()); fate["topabsP"] += int((~refl & ~trans & ~a_).sum())
                alive[ht[~refl]] = False
            hb = bz[~up]
            if hb.size:
                z[hb] = 0.0
                bbin = np.minimum(np.searchsorted(edges, np.arccos(np.clip(-uz[hb], 0, 1))) - 1, nth - 1)
                a_ = (np.abs(x[hb]) <= h) & (np.abs(y[hb]) <= h)
                Rv = np.where(a_, tA["Rb"][pol[hb], il[hb], bbin], tP["Rb"][pol[hb], il[hb], bbin])
                ok = rng.random(hb.size) < Rv
                uz[hb[ok]] *= -1
                fate["bot"] += int((~ok).sum()); alive[hb[~ok]] = False
        sc = ids[kind == 3]
        if sc.size:
            ct = rc._hg_cos(g, rng.random(sc.size)); st = np.sqrt(np.clip(1 - ct * ct, 0, 1))
            ph = 2 * np.pi * rng.random(sc.size); cp, sp = np.cos(ph), np.sin(ph)
            ux0, uy0, uz0 = ux[sc], uy[sc], uz[sc]
            tmp = np.sqrt(np.clip(1 - uz0 * uz0, 1e-12, 1)); big = np.abs(uz0) > 0.99999
            nx = np.where(big, st * cp, st * (ux0 * uz0 * cp - uy0 * sp) / tmp + ux0 * ct)
            ny = np.where(big, st * sp, st * (uy0 * uz0 * cp + ux0 * sp) / tmp + uy0 * ct)
            nz = np.where(big, np.sign(uz0) * ct, -st * cp * tmp + uz0 * ct)
            nr = np.sqrt(nx * nx + ny * ny + nz * nz)
            ux[sc], uy[sc], uz[sc] = nx / nr, ny / nr, nz / nr
            pol[sc] = (rng.random(sc.size) < 0.5).astype(pol.dtype)
    fate["lost"] = int(alive.sum())
    sc_ = total_in / n_photons
    return dict(arrA=arrA * sc_, arrP=arrP * sc_, fate=fate, scale=sc_, n=n_photons)


if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "prop3_narrow528_ext"
    d = an.load(tag)
    pm, pf, _ = opt.preset_mats(d.get("preset"))
    ctx, s0 = opt.make_ctx(d["emitter"], mats=pm, step=4.0)
    P = an._restore(d["params"])
    run = run_proposed(ctx, P, s0, n_th=91, du=1e-3)
    th = run["th"]; nm = run["P"]["n_m"]
    mir = ctx.nk(ctx.mats["mirror"])
    keep = len(sys.argv) > 3 and sys.argv[3] == "keep_top"
    tA = run["tables"]; tP = apron_tables(ctx, run["P"], th, mir, keep_ag=keep)
    T0A = 0.5 * (tA["Tt"][0, :, 0] + tA["Tt"][1, :, 0]); T0P = 0.5 * (tP["Tt"][0, :, 0] + tP["Tt"][1, :, 0])
    cases = [  # (label, d_um, musA, musP, g, L, a)
        ("diffusive d2 a15", 2.0, 2.5, 2.5, 0.75, 50, 15),
        ("clear-core d2 a15", 2.0, 0.02, 2.5, 0.75, 50, 15),
        ("diffusive d5 a15", 5.0, 1.0, 1.0, 0.75, 50, 15),
        ("clear-core d5 a15", 5.0, 0.02, 1.0, 0.75, 50, 15),
        ("clear-core d5 a10 L50", 5.0, 0.02, 1.0, 0.75, 50, 10),
        ("diffusive d2 a5 L17", 2.0, 2.5, 2.5, 0.75, 17, 5),
        ("clear-core d2 a5 L17", 2.0, 0.02, 2.5, 0.75, 17, 5),
        ("clear-core d5 full-cell(a=L)", 5.0, 0.02, 1.0, 0.75, 50, 50),
    ]
    only = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    rows = []
    for i, (lab, dd, mA, mP, g, L, a) in enumerate(cases):
        if only and str(i) not in only:
            continue
        mc = mc_cell2(tA, tP, th, ctx.wl, run["Pin"], run["qw"], dd, mA, mP, g, L, a)
        ce = ce_from_arrivals(mc["arrA"], T0A, ctx.wl, ctx.V, nm, th) + ce_from_arrivals(mc["arrP"], T0P, ctx.wl, ctx.V, nm, th)
        f = mc["fate"]; n = mc["n"]
        row = dict(label=lab, d=dd, musA=mA, musP=mP, L=L, a=a, fill=(a / L) ** 2, CE_total=float(ce + run["summary"]["CE_dir"]),
                   eqe_total=float((f["escA"] + f["escP"]) * mc["scale"] + run["summary"]["eqe_dir"]),
                   frac_from_apron=f["escP"] / max(1, f["escA"] + f["escP"]), loss_topA=f["topabsA"] / n,
                   loss_topP=f["topabsP"] / n, loss_bot=f["bot"] / n, loss_wall=f["wall"] / n, lost=f["lost"] / n)
        rows.append(row)
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    print("analytic laterally-infinite:", round(run["summary"]["CE"], 1), round(run["summary"]["eqe"], 3))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", f"apron2_{tag}{'_keeptop' if keep else ''}.json"), "w") as fh:
        json.dump(dict(tag=tag, apron_keeps_top_electrode=keep, rows=rows, base_CE=run["summary"]["CE"]), fh, indent=1)
