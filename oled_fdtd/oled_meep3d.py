"""3D version of the minimal OLED power-budget FDTD (see oled_meep.py).

Same stack, now a true 3D cell so results are directly comparable with the
classical CPS / transfer-matrix dipole model (point dipole physics).

Axes: y = vertical (stack normal), x/z = lateral.
Dipole orientations: horizontal = Ex (Ez equivalent by symmetry, weighted x2
in the isotropic average), vertical = Ey.  Isotropic avg = (2*Ex + Ey)/3.
"""

import argparse
import json

import meep as mp
from meep.materials import Al

UM = 1.0
LAM = 0.550 * UM
FCEN = 1 / LAM
DF = 0.2 * FCEN

N_GLASS = 1.5
N_ITO = 1.9
N_ORG = 1.75

T_ITO = 0.100 * UM
T_ORG = 0.100 * UM
T_AL = 0.200 * UM
T_GLASS = 1.0 * UM

DPML = 0.4 * UM
SX = 3.0 * UM            # lateral extent (x and z)


def run(pol, resolution=50, courant=0.5):
    sy = DPML + T_GLASS + T_ITO + T_ORG + T_AL + DPML
    cell = mp.Vector3(SX + 2 * DPML, sy, SX + 2 * DPML)

    y0 = -sy / 2
    y_al0 = y0 + DPML
    y_al1 = y_al0 + T_AL
    y_org1 = y_al1 + T_ORG
    y_ito1 = y_org1 + T_ITO

    def slab(y_lo, y_hi, mat):
        return mp.Block(center=mp.Vector3(0, (y_lo + y_hi) / 2, 0),
                        size=mp.Vector3(mp.inf, y_hi - y_lo, mp.inf),
                        material=mat)

    geometry = [
        slab(y_ito1, -y0, mp.Medium(index=N_GLASS)),
        slab(y_org1, y_ito1, mp.Medium(index=N_ITO)),
        slab(y_al1, y_org1, mp.Medium(index=N_ORG)),
        slab(y_al0, y_al1, Al),
    ]

    y_src = (y_al1 + y_org1) / 2
    comp = {"Ex": mp.Ex, "Ey": mp.Ey, "Ez": mp.Ez}[pol]
    src = [mp.Source(mp.GaussianSource(FCEN, fwidth=DF),
                     component=comp, center=mp.Vector3(0, y_src, 0))]

    # vertical PML; lateral Absorber (SPP-in-PML instability, see 2D version)
    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        sources=src,
                        boundary_layers=[mp.PML(DPML, direction=mp.Y),
                                         mp.Absorber(DPML, direction=mp.X),
                                         mp.Absorber(DPML, direction=mp.Z)],
                        resolution=resolution,
                        Courant=courant)

    # ---- flux surfaces -------------------------------------------------
    d = 0.06 * UM
    def fr(center, size, direction, weight):
        return sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
            center=center, size=size, direction=direction, weight=weight))

    c = mp.Vector3(0, y_src, 0)
    box = [
        fr(c + mp.Vector3(0, d / 2, 0), mp.Vector3(d, 0, d), mp.Y, +1),
        fr(c - mp.Vector3(0, d / 2, 0), mp.Vector3(d, 0, d), mp.Y, -1),
        fr(c + mp.Vector3(d / 2, 0, 0), mp.Vector3(0, d, d), mp.X, +1),
        fr(c - mp.Vector3(d / 2, 0, 0), mp.Vector3(0, d, d), mp.X, -1),
        fr(c + mp.Vector3(0, 0, d / 2), mp.Vector3(d, d, 0), mp.Z, +1),
        fr(c - mp.Vector3(0, 0, d / 2), mp.Vector3(d, d, 0), mp.Z, -1),
    ]

    glass_plane = fr(mp.Vector3(0, y_ito1 + 0.15 * UM, 0),
                     mp.Vector3(SX, 0, SX), mp.Y, +1)

    y_wg_lo = y_al1 - 0.05 * UM
    y_wg_hi = y_ito1 + 0.10 * UM
    ywg = (y_wg_lo + y_wg_hi) / 2
    hwg = y_wg_hi - y_wg_lo
    xedge = SX / 2 - 0.05 * UM
    wg = [
        fr(mp.Vector3(-xedge, ywg, 0), mp.Vector3(0, hwg, SX), mp.X, -1),
        fr(mp.Vector3(+xedge, ywg, 0), mp.Vector3(0, hwg, SX), mp.X, +1),
        fr(mp.Vector3(0, ywg, -xedge), mp.Vector3(SX, hwg, 0), mp.Z, -1),
        fr(mp.Vector3(0, ywg, +xedge), mp.Vector3(SX, hwg, 0), mp.Z, +1),
    ]

    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        5, comp, mp.Vector3(0, y_src, 0), 1e-6))

    p_tot = sum(mp.get_fluxes(f)[0] for f in box)
    p_glass = mp.get_fluxes(glass_plane)[0]
    p_wg = sum(mp.get_fluxes(f)[0] for f in wg)

    out = dict(pol=pol, P_total=p_tot,
               frac_glass=p_glass / p_tot,
               frac_wg=p_wg / p_tot)
    out["frac_metal"] = 1 - out["frac_glass"] - out["frac_wg"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=int, default=50)
    args = ap.parse_args()

    results = {}
    for pol in ("Ex", "Ey"):        # Ez == Ex by symmetry
        r = run(pol, resolution=args.res)
        results[pol] = r
        print("[3D %s res=%d]  glass %.3f  wg %.3f  metal %.3f"
              % (pol, args.res, r["frac_glass"], r["frac_wg"], r["frac_metal"]))

    iso = {k: (2 * results["Ex"][k] + results["Ey"][k]) / 3
           for k in ("frac_glass", "frac_wg", "frac_metal")}
    print("[3D isotropic (2*Ex+Ey)/3]  glass %.3f  wg %.3f  metal %.3f"
          % (iso["frac_glass"], iso["frac_wg"], iso["frac_metal"]))

    with open("oled_meep3d_res%d.json" % args.res, "w") as f:
        json.dump({"per_pol": results, "isotropic": iso,
                   "lambda_um": LAM, "resolution": args.res}, f, indent=2)
    print("OLED3D-DONE")


if __name__ == "__main__":
    main()
