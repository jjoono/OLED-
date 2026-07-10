"""Minimal 2D FDTD simulation of a bottom-emitting OLED with Meep (open source).

Goal: power budget of dipole emission in the simplest credible OLED stack,
at a single wavelength (lambda = 550 nm).

Stack (top of cell -> bottom), units of 1 um:

        PML
        glass substrate   n = 1.5   (semi-infinite -> ends in PML)
        ITO               n = 1.9   100 nm   (lossless, simplest)
        organic (EML)     n = 1.75  100 nm   <- point dipole in the middle
        Al cathode        Rakic Drude-Lorentz (meep.materials.Al)
        PML

Outputs, per dipole orientation (Ex / Ey / Ez) and isotropic average:
    frac_glass  : power escaping upward into the glass substrate
                  ("substrate-coupled" light; the part that can reach air)
    frac_wg     : power flowing sideways in the ITO/organic waveguide
    frac_metal  : remainder = absorption in the Al (incl. SPP)   (1 - above)

Method: narrowband Gaussian pulse around f0 = 1/0.55; run until fields decay;
Poynting flux through planes, normalised by the total power emitted by the
dipole (flux through a small box enclosing the source).
"""

import argparse
import json

import meep as mp
from meep.materials import Al

# ---------------------------------------------------------------- geometry
UM = 1.0
LAM = 0.550 * UM          # vacuum wavelength
FCEN = 1 / LAM
DF = 0.2 * FCEN           # narrowband pulse

N_GLASS = 1.5
N_ITO = 1.9
N_ORG = 1.75

T_ITO = 0.100 * UM
T_ORG = 0.100 * UM
T_AL = 0.200 * UM         # opaque at 550 nm
T_GLASS = 1.2 * UM        # "semi-infinite": thick slab ending in PML

DPML = 0.5 * UM
SX = 6.0 * UM             # lateral size (waveguide power measured near edges)

RES = 100                 # pixels / um  -> 10 nm grid


def build(pol, resolution=RES, courant=0.5):
    sy = DPML + T_GLASS + T_ITO + T_ORG + T_AL + DPML
    cell = mp.Vector3(SX + 2 * DPML, sy)

    # y coordinates (cell centred at 0): put organic/ITO interface near centre
    y0 = -sy / 2
    y_al0 = y0 + DPML
    y_al1 = y_al0 + T_AL
    y_org1 = y_al1 + T_ORG
    y_ito1 = y_org1 + T_ITO
    # glass fills the rest up to the top PML

    def slab(y_lo, y_hi, mat):
        return mp.Block(center=mp.Vector3(0, (y_lo + y_hi) / 2),
                        size=mp.Vector3(mp.inf, y_hi - y_lo, mp.inf),
                        material=mat)

    geometry = [
        slab(y_ito1, -y0, mp.Medium(index=N_GLASS)),
        slab(y_org1, y_ito1, mp.Medium(index=N_ITO)),
        slab(y_al1, y_org1, mp.Medium(index=N_ORG)),
        slab(y_al0, y_al1, Al),
    ]

    y_src = (y_al1 + y_org1) / 2          # middle of the organic layer
    comp = {"Ex": mp.Ex, "Ey": mp.Ey, "Ez": mp.Ez}[pol]
    src = [mp.Source(mp.GaussianSource(FCEN, fwidth=DF),
                     component=comp, center=mp.Vector3(0, y_src))]

    sim = mp.Simulation(cell_size=cell,
                        geometry=geometry,
                        sources=src,
                        boundary_layers=[mp.PML(DPML)],
                        resolution=resolution,
                        Courant=courant,
                        force_complex_fields=False)

    return sim, dict(y_src=y_src, y_ito1=y_ito1, y_al1=y_al1, y_org1=y_org1,
                     sy=sy, y0=y0)


def run(pol, resolution=RES, quiet=True, courant=0.5):
    sim, g = build(pol, resolution, courant)
    y_src, y_ito1 = g["y_src"], g["y_ito1"]

    # total emitted power: small box around the dipole (inside the organic)
    d = 0.06 * UM
    box = [
        sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
            center=mp.Vector3(0, y_src + d / 2), size=mp.Vector3(d, 0),
            direction=mp.Y, weight=+1)),
        sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
            center=mp.Vector3(0, y_src - d / 2), size=mp.Vector3(d, 0),
            direction=mp.Y, weight=-1)),
        sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
            center=mp.Vector3(+d / 2, y_src), size=mp.Vector3(0, d),
            direction=mp.X, weight=+1)),
        sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
            center=mp.Vector3(-d / 2, y_src), size=mp.Vector3(0, d),
            direction=mp.X, weight=-1)),
    ]

    # power into the glass: plane just above the ITO/glass interface
    glass_plane = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
        center=mp.Vector3(0, y_ito1 + 0.15 * UM),
        size=mp.Vector3(SX, 0), direction=mp.Y, weight=+1))

    # waveguided power: vertical planes near the lateral edges, spanning
    # organic + ITO (+ a little of the adjacent media to catch the mode tails)
    y_wg_lo = g["y_al1"] - 0.05 * UM
    y_wg_hi = y_ito1 + 0.10 * UM
    wg_l = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
        center=mp.Vector3(-SX / 2 + 0.1 * UM, (y_wg_lo + y_wg_hi) / 2),
        size=mp.Vector3(0, y_wg_hi - y_wg_lo), direction=mp.X, weight=-1))
    wg_r = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
        center=mp.Vector3(+SX / 2 - 0.1 * UM, (y_wg_lo + y_wg_hi) / 2),
        size=mp.Vector3(0, y_wg_hi - y_wg_lo), direction=mp.X, weight=+1))

    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        5, {"Ex": mp.Ex, "Ey": mp.Ey, "Ez": mp.Ez}[pol],
        mp.Vector3(0, y_src), 1e-6))

    p_tot = sum(mp.get_fluxes(f)[0] for f in box)
    p_glass = mp.get_fluxes(glass_plane)[0]
    p_wg = mp.get_fluxes(wg_l)[0] + mp.get_fluxes(wg_r)[0]

    out = dict(pol=pol,
               P_total=p_tot,
               frac_glass=p_glass / p_tot,
               frac_wg=p_wg / p_tot)
    out["frac_metal"] = 1 - out["frac_glass"] - out["frac_wg"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=int, default=RES)
    args = ap.parse_args()

    results = {}
    for pol in ("Ez", "Ex", "Ey"):
        r = run(pol, resolution=args.res)
        results[pol] = r
        print("[%s]  glass %.3f   waveguide %.3f   metal(SPP+abs) %.3f"
              % (pol, r["frac_glass"], r["frac_wg"], r["frac_metal"]))

    iso = {k: sum(results[p][k] for p in results) / 3
           for k in ("frac_glass", "frac_wg", "frac_metal")}
    print("[isotropic 1/3 avg]  glass %.3f   waveguide %.3f   metal %.3f"
          % (iso["frac_glass"], iso["frac_wg"], iso["frac_metal"]))

    with open("oled_meep_result.json", "w") as f:
        json.dump({"per_pol": results, "isotropic": iso,
                   "lambda_um": LAM, "resolution": args.res}, f, indent=2)


if __name__ == "__main__":
    main()
