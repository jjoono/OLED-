"""
optimize.py -- global optimisation (differential evolution) of on-axis current efficiency
for the proposed angular-recycling TEOLED and for a conventional strong-microcavity TEOLED
with identical organic stack, emitter and top anode.

usage: python3 optimize.py <job> [emitter] [extra json overrides]
jobs : ref | prop | propfree
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.optimize import differential_evolution

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oledsim import MaterialDB  # noqa: E402
from devices import Ctx, run_proposed, run_reference  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DB = MaterialDB(os.path.join(HERE, "..", "nk_JH_total.mat"))

EMITTERS = {
    "narrow528": dict(kind="narrow", peak=528, fwhm=25, shoulder=0.15),
    "narrow532_18": dict(kind="narrow", peak=532, fwhm=18, shoulder=0.10),
    "irppy2acac": dict(kind="l_I_Irppy2acac"),
    "vdabna_hf": dict(kind="l_mCPCN_SpiroAC_20wt_v_DABNA_1_5wt"),   # measured blue hyperfluorescence, 475 nm / 23 nm
    "red_narrow": dict(kind="narrow", peak=615, fwhm=30, shoulder=0.15),
    "narrow532_15": dict(kind="narrow", peak=532, fwhm=15, shoulder=0.0),   # ultra-narrow, no vibronic shoulder
    "blue_narrow": dict(kind="narrow", peak=470, fwhm=20, shoulder=0.10),       # nu-DABNA-like EL (469 nm / 18 nm)
    "narrow528_40": dict(kind="narrow", peak=528, fwhm=40, shoulder=0.15),
}


def preset_mats(preset):
    """material / fixed-parameter presets for A' engineering studies."""
    mats, fixed, ag_b = {}, {}, None
    for tok in (preset or "").split("+"):
        if not tok:
            continue
        if tok == "bulkAg":        # thin Ag with bulk-like optics (seeded / alloyed growth limit)
            mats["ag_top"] = "Ag_McPeak"
        elif tok == "Ag12":        # measured 12-nm Ag (n~0.09 @550)
            mats["ag_top"] = "l_Ag_12nm_test5"
        elif tok == "lowkITO":     # ITO with extinction coefficient scaled to 30 %
            ito = DB.nk("l_ITO", np.arange(400.0, 801.0))
            mats["tco"] = ("lowkITO", ito.real + 0.3j * ito.imag)
        elif tok == "IZO":
            mats["tco"] = "l_IZO_SJ_RTP"
        elif tok == "ODR":         # 800-nm SiO2 spacer between SL and the Ag mirror (>=500 nm avoids Otto-SPP coupling)
            fixed["spacer"] = ("l_SiO2_SJ_RTP", 800.0)
        elif tok == "lowHTL":      # hypothetical low-index HTL (n = 1.55)
            mats["htl"] = (1.55, 1.55)
        elif tok == "IZOtop":      # sputtered IZO top anode instead of thin Ag (damage-free limit)
            mats["ag_top"] = "l_IZO_SJ_RTP"
            ag_b = (20.0, 150.0)
        else:
            raise ValueError(tok)
    return mats, fixed, ag_b


def make_ctx(emitter, step=4.0, mats=None, full=False):
    wl_all = np.arange(400.0, 701.0, 2.0)
    c0 = Ctx(DB, wl=wl_all, mats=mats)
    s = c0.emitter(**EMITTERS[emitter])
    if full:
        return c0, s
    keep = s > 2e-3 * s.max()
    lo, hi = wl_all[keep].min(), wl_all[keep].max()
    wl = np.arange(lo, hi + 0.1, step)
    c = Ctx(DB, wl=wl, mats=mats)
    return c, c.emitter(**EMITTERS[emitter])


# ----------------------------------------------------------------------------
G = {}


def _ref_params(x):
    d_etl, d_htl, d_ag, d_cap = x
    return dict(d_etl=d_etl, d_htl=d_htl, d_ag=d_ag, cap=("cpl", d_cap))


def _prop_params(x, N):
    d_etl, d_htl, d_ag, dH1, dH, dL, n_m = x
    return dict(d_etl=d_etl, d_htl=d_htl, d_ag=d_ag, n_m=n_m, dbr=dict(N=N, dH=dH, dL=dL, dH1=dH1))


def _propfree_params(x, N):
    d_etl, d_htl, d_ag, n_m = x[:4]
    return dict(d_etl=d_etl, d_htl=d_htl, d_ag=d_ag, n_m=n_m, dbr=list(x[4:4 + 2 * N + 1]))


def _refdbr_params(x, N):
    d_etl, d_htl, d_ag, dH1, dH, dL = x
    return dict(d_etl=d_etl, d_htl=d_htl, d_ag=d_ag, cap=None, dbr=dict(N=N, dH=dH, dL=dL, dH1=dH1))


def obj_refdbr(x):
    r = run_reference(G["ctx"], {**G["fixed"], **_refdbr_params(x, G["N"])}, G["s0"], du=1e-3)
    return -r["summary"]["CE"]


def _cie_y(ctx, spec_ph):
    xb, yb, zb = ctx.db.xyz_bar(ctx.wl)
    e = spec_ph * (1239.841984 / ctx.wl)
    X, Y, Z = np.sum(xb * e), np.sum(yb * e), np.sum(zb * e)
    return Y / (X + Y + Z)


def _score(r, kind):
    """objective: on-axis CE, or blue index BI = CE / CIE_y of the on-axis spectrum (G['objective'] == 'BI')."""
    ce = r["summary"]["CE"]
    if G.get("objective", "CE") != "BI":
        return -ce
    if "out" in r:  # proposed
        sp = r["out"]["L_dir0"] + r["out"]["L_rec0"]
    else:
        from oledsim import radiant_intensity
        Is, Ip = radiant_intensity(r["res"], "top", 1.0, np.array([0.0]))
        sp = r["qw"] * (Is + Ip)[:, 0]
    return -ce / max(_cie_y(G["ctx"], sp), 1e-3)


def obj_ref(x):
    r = run_reference(G["ctx"], {**G["fixed"], **_ref_params(x)}, G["s0"], du=1e-3)
    return _score(r, "ref")


def _prophyb_params(x, N, Nb):
    d_etl, d_htl, d_ag, dH1, dH, dL, n_m, d_tco, bH, bL = x
    return dict(d_etl=d_etl, d_htl=d_htl, d_ag=d_ag, n_m=n_m, d_tco=d_tco,
                dbr=dict(N=N, dH=dH, dL=dL, dH1=dH1), bdbr=dict(N=Nb, dH=bH, dL=bL, dH1=bH))


def obj_prophyb(x):
    r = run_proposed(G["ctx"], {**G["fixed"], **_prophyb_params(x, G["N"], G["Nb"])}, G["s0"], n_th=91, du=1e-3,
                     a_sl_pass=G.get("a_sl", 0.0))
    return -r["summary"]["CE"]


def obj_prop(x):
    r = run_proposed(G["ctx"], {**G["fixed"], **_prop_params(x, G["N"])}, G["s0"], n_th=91, du=1e-3,
                     a_sl_pass=G.get("a_sl", 0.0))
    return _score(r, "prop")


def obj_propfree(x):
    r = run_proposed(G["ctx"], {**G["fixed"], **_propfree_params(x, G["N"])}, G["s0"], n_th=91, du=1e-3,
                     a_sl_pass=G.get("a_sl", 0.0))
    return -r["summary"]["CE"]


def run(job, emitter="narrow528", N=3, fixed=None, mats=None, maxiter=120, popsize=15, seed=0,
        a_sl=0.0, d_ag_min=10.0, n_m_bounds=(1.6, 2.0), preset=None, Nb=2, step=4.0, objective="CE"):
    pm, pf, ag_b = preset_mats(preset)
    mats = {**pm, **(mats or {})}
    fixed = {**pf, **(fixed or {})}
    ctx, s0 = make_ctx(emitter, mats=mats, step=step)
    G.update(ctx=ctx, s0=s0, fixed=fixed or {}, N=N, a_sl=a_sl, Nb=Nb, objective=objective)
    t0 = time.time()
    if job == "ref":
        bounds = [(10, 400), (20, 400), (d_ag_min, 40), (0, 160)]
        f = obj_ref
    elif job == "refdbr":
        bounds = [(10, 320), (20, 320), (d_ag_min, 30), (0, 220), (20, 220), (20, 260)]
        f = obj_refdbr
    elif job == "prop":
        bounds = [(10, 400), (20, 500), ag_b or (d_ag_min, 25), (0, 220), (20, 220), (20, 300), n_m_bounds]
        f = obj_prop
    elif job == "prophyb":
        bounds = [(10, 320), (20, 320), ag_b or (d_ag_min, 25), (0, 220), (20, 220), (20, 260), n_m_bounds,
                  (20, 150), (20, 160), (40, 260)]
        f = obj_prophyb
    else:
        bounds = [(10, 320), (20, 320), ag_b or (d_ag_min, 25), n_m_bounds] + [(0, 260)] * (2 * N + 1)
        f = obj_propfree
    res = differential_evolution(f, bounds, maxiter=maxiter, popsize=popsize, seed=seed, tol=1e-6,
                                 mutation=(0.5, 1.0), recombination=0.8, polish=False, workers=4,
                                 updating="deferred")
    # final evaluation on a fine grid
    ctxf, s0f = make_ctx(emitter, mats=mats, full=True)
    if job == "ref":
        P = {**(fixed or {}), **_ref_params(res.x)}
        fin = run_reference(ctxf, P, s0f, du=2.5e-4)["summary"]
    elif job == "refdbr":
        P = {**(fixed or {}), **_refdbr_params(res.x, N)}
        fin = run_reference(ctxf, P, s0f, du=2.5e-4)["summary"]
    elif job == "prop":
        P = {**(fixed or {}), **_prop_params(res.x, N)}
        fin = run_proposed(ctxf, P, s0f, du=2.5e-4, a_sl_pass=a_sl)["summary"]
    elif job == "prophyb":
        P = {**(fixed or {}), **_prophyb_params(res.x, N, Nb)}
        fin = run_proposed(ctxf, P, s0f, du=2.5e-4, a_sl_pass=a_sl)["summary"]
    else:
        P = {**(fixed or {}), **_propfree_params(res.x, N)}
        fin = run_proposed(ctxf, P, s0f, du=2.5e-4, a_sl_pass=a_sl)["summary"]
    out = dict(job=job, emitter=emitter, N=N, x=list(map(float, res.x)), fun=float(res.fun), objective=objective,
               params=_jsonable(P),
               final=fin, fixed=_jsonable(fixed or {}), preset=preset, a_sl=a_sl, nfev=int(res.nfev),
               seconds=time.time() - t0)
    return out


def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


if __name__ == "__main__":
    job = sys.argv[1]
    emitter = sys.argv[2] if len(sys.argv) > 2 else "narrow528"
    kw = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    tag = kw.pop("tag", f"{job}_{emitter}")
    out = run(job, emitter, **kw)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", f"opt_{tag}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
