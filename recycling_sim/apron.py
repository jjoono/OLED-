"""
apron.py -- spatial + angular recycling: a pixel cell (square, side L) whose central square
(side a) is the electrically active OLED, surrounded by an electrode-free, low-loss "apron"
(same SL + mirror below, same DBR above, but no top thin Ag and no bottom TCO), bounded by
reflective banks (reflectance R_wall).  Light emitted downward in the active area spreads
laterally in the SL and is recycled; in the apron the round-trip loss is much smaller.

Figure of merit: on-axis luminous intensity per current (cd/A in the area-independent sense,
i.e. total on-axis intensity of the pixel / pixel current).  Panel luminance at fixed pixel pitch
scales the same way.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import recycling as rc  # noqa: E402
from oledsim import KM, HC_E, Medium, stack_RTA  # noqa: E402
from devices import run_proposed, proposed_layers  # noqa: E402
import optimize as opt  # noqa: E402
import analyze as an  # noqa: E402


def mc_cell(tA, tP, th_deg, wl, Pin, qw, n_m, d_um, mus, g, L_um, a_um, R_wall=0.95,
            n_photons=150000, max_steps=40000, seed=3, shifts=None):
    """2-zone photon MC.  tA/tP: boundary tables (Rt, Tt, Rb) for active / apron zones.
    Returns per-exciton weighted near-normal arrivals (for the radiance estimator) and fates."""
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
    x = (rng.random(n_photons) - 0.5) * a_um
    y = (rng.random(n_photons) - 0.5) * a_um
    z = np.full(n_photons, d_um)
    alive = np.ones(n_photons, bool)
    mut = mus
    edges = np.deg2rad(np.linspace(0, 90, nth + 1))
    arrA = np.zeros((2, nwl, nth)); arrP = np.zeros((2, nwl, nth))
    fate = dict(escA=0, escP=0, topabsA=0, topabsP=0, bot=0, wall=0, lost=0)
    half = L_um / 2
    for step in range(max_steps):
        ids = np.nonzero(alive)[0]
        if ids.size == 0:
            break
        uzi = uz[ids]
        dist_z = np.where(uzi > 0, (d_um - z[ids]) / np.maximum(uzi, 1e-300), z[ids] / np.maximum(-uzi, 1e-300))
        # distance to lateral walls
        with np.errstate(divide="ignore", invalid="ignore"):
            dx = np.where(ux[ids] > 0, (half - x[ids]) / ux[ids], np.where(ux[ids] < 0, (-half - x[ids]) / ux[ids], np.inf))
            dy = np.where(uy[ids] > 0, (half - y[ids]) / uy[ids], np.where(uy[ids] < 0, (-half - y[ids]) / uy[ids], np.inf))
        s = -np.log(rng.random(ids.size)) / mut
        dmin = np.minimum(np.minimum(dist_z, dx), np.minimum(dy, s))
        x[ids] += ux[ids] * dmin; y[ids] += uy[ids] * dmin; z[ids] += uz[ids] * dmin
        # walls
        hw = (dmin == dx) | (dmin == dy)
        wids = ids[hw]
        if wids.size:
            ok = rng.random(wids.size) < R_wall
            fate["wall"] += int((~ok).sum()); alive[wids[~ok]] = False
            w2 = wids[ok]
            hx = np.abs(np.abs(x[w2]) - half) < 1e-9
            ux[w2[hx]] *= -1; uy[w2[~hx]] *= -1
            x[w2] = np.clip(x[w2], -half + 1e-9, half - 1e-9); y[w2] = np.clip(y[w2], -half + 1e-9, half - 1e-9)
        rest = ids[~hw]; dm = dmin[~hw]; dz = dist_z[~hw]; ss = s[~hw]
        hz = (dm == dz)
        bz = rest[hz]
        if bz.size:
            up = uz[bz] > 0
            ht = bz[up]
            if ht.size:
                z[ht] = d_um
                tb = np.minimum(np.searchsorted(edges, np.arccos(np.clip(uz[ht], 0, 1))) - 1, nth - 1)
                act = (np.abs(x[ht]) <= a_um / 2) & (np.abs(y[ht]) <= a_um / 2)
                Rv = np.where(act, tA["Rt"][pol[ht], il[ht], tb], tP["Rt"][pol[ht], il[ht], tb])
                Tv = np.where(act, tA["Tt"][pol[ht], il[ht], tb], tP["Tt"][pol[ht], il[ht], tb])
                np.add.at(arrA, (pol[ht[act]], il[ht[act]], tb[act]), 1.0)
                np.add.at(arrP, (pol[ht[~act]], il[ht[~act]], tb[~act]), 1.0)
                xi = rng.random(ht.size)
                refl = xi < Rv; trans = (~refl) & (xi < Rv + Tv)
                uz[ht[refl]] *= -1
                fate["escA"] += int((trans & act).sum()); fate["escP"] += int((trans & ~act).sum())
                fate["topabsA"] += int((~refl & ~trans & act).sum()); fate["topabsP"] += int((~refl & ~trans & ~act).sum())
                alive[ht[~refl]] = False
            hb = bz[~up]
            if hb.size:
                z[hb] = 0.0
                bb = np.minimum(np.searchsorted(edges, np.arccos(np.clip(-uz[hb], 0, 1))) - 1, nth - 1)
                act = (np.abs(x[hb]) <= a_um / 2) & (np.abs(y[hb]) <= a_um / 2)
                Rv = np.where(act, tA["Rb"][pol[hb], il[hb], bb], tP["Rb"][pol[hb], il[hb], bb])
                ok = rng.random(hb.size) < Rv
                uz[hb[ok]] *= -1
                fate["bot"] += int((~ok).sum()); alive[hb[~ok]] = False
        sc = rest[~hz]
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


def ce_from_arrivals(arr, T0, wl, V, n_m, th_deg, cone_deg=10.0):
    cm = np.arcsin(np.sin(np.deg2rad(cone_deg)) / n_m)
    edges = np.deg2rad(np.linspace(0, 90, len(th_deg) + 1))
    nb = max(1, np.searchsorted(edges, cm) - 1)
    Lm = arr[:, :, :nb].sum(axis=(0, 2)) / (np.pi * np.sin(edges[nb]) ** 2)
    return KM * np.sum(V * (HC_E / wl) * T0 * Lm / n_m ** 2)


def apron_tables(ctx, P, th, mirror_nk, keep_ag=False):
    _, _, _, top, bot = proposed_layers(ctx, P)
    drop = {"tco", "eil"} | (set() if keep_ag else {"ag_top"})
    topP = [l for l in top if l.name not in drop]
    return rc.boundary_tables(topP, P["n_m"] * ctx.ones, ctx.wl, th, bot_layers=bot, mirror_nk=mirror_nk)


if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "prop3_narrow528_ext"
    d = an.load(tag)
    pm, pf, _ = opt.preset_mats(d.get("preset"))
    ctx, s0 = opt.make_ctx(d["emitter"], mats=pm, step=4.0)
    P = an._restore(d["params"])
    run = run_proposed(ctx, P, s0, n_th=91, du=1e-3)
    th = run["th"]; nm = P["n_m"] if "n_m" in P else 1.85
    run["P"]["n_m"] = run["P"].get("n_m", nm)
    mir = ctx.nk(ctx.mats["mirror"])
    tA = run["tables"]
    rows = []
    for keep_ag in (False, True):
        tP = apron_tables(ctx, run["P"], th, mir, keep_ag=keep_ag)
        T0A = 0.5 * (tA["Tt"][0, :, 0] + tA["Tt"][1, :, 0]); T0P = 0.5 * (tP["Tt"][0, :, 0] + tP["Tt"][1, :, 0])
        w = 0.5 * run["wL"][None, None, :]
        pP = float(np.sum(tP["Tt"] * w, axis=(0, 2)).mean()); aP = float(np.sum(tP["At"] * w, axis=(0, 2)).mean())
        for (L, a, Rw) in [(50, 50, 0.95), (50, 35, 0.95), (50, 25, 0.95), (50, 15, 0.95), (50, 25, 0.99),
                           (100, 25, 0.95), (30, 15, 0.95)]:
            mc = mc_cell(tA, tP, th, ctx.wl, run["Pin"], run["qw"], run["P"]["n_m"], 2.0, 2.5, 0.75, L, a, R_wall=Rw,
                         n_photons=120000)
            ce = ce_from_arrivals(mc["arrA"], T0A, ctx.wl, ctx.V, run["P"]["n_m"], th) + \
                 ce_from_arrivals(mc["arrP"], T0P, ctx.wl, ctx.V, run["P"]["n_m"], th)
            f = mc["fate"]; n = mc["n"]
            eqe_rec = (f["escA"] + f["escP"]) * mc["scale"]
            row = dict(apron_has_Ag=keep_ag, L=L, a=a, fill=(a / L) ** 2, R_wall=Rw, CE_rec=float(ce),
                       CE_total=float(ce + run["summary"]["CE_dir"]), eqe_total=float(eqe_rec + run["summary"]["eqe_dir"]),
                       frac_from_apron=f["escP"] / max(1, f["escA"] + f["escP"]),
                       loss_topA=f["topabsA"] / n, loss_topP=f["topabsP"] / n, loss_bot=f["bot"] / n, loss_wall=f["wall"] / n,
                       apron_p=pP, apron_a_top=aP, apron_T0peak=float(T0P[np.argmax(s0)]))
            rows.append(row)
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
    print("reference analytic (laterally infinite active):", round(run["summary"]["CE"], 1), round(run["summary"]["eqe"], 3))
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", f"apron_{tag}.json"), "w") as fh:
        json.dump(dict(tag=tag, rows=rows, base_CE=run["summary"]["CE"], base_eqe=run["summary"]["eqe"]), fh, indent=1)
