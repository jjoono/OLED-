"""Substrate -> device plane-wave R/T for a periodically patterned (island) layer.

2D FDTD (Meep). A plane wave is launched in the glass substrate toward the
device; we measure the reflected power back into the glass and the transmitted
power into the superstrate. When the patterned layer is filled to 100 %
(a continuous slab) the result must reproduce the analytic TMM R/T of a
planar glass / film / air stack -- that is the validation. With islands
(fill < 1) the layer diffracts and TMM no longer applies; FDTD still gives the
correct total R/T (and R + T = 1 with lossless media, an energy check).

Geometry (y = propagation, x = periodic with period Lambda):

    PML
    air  (superstrate, n=1.0)          <- transmission monitor
    patterned layer, thickness d:
        island  n=1.75  width = fill*Lambda
        gap     n=1.00  (air)
    glass (substrate, n=1.5)           <- reflection monitor, plane-wave source
    PML

Single wavelength lambda = 550 nm, normal incidence, Ez polarization.
"""

import argparse
import numpy as np
import meep as mp

LAM = 0.550
FCEN = 1 / LAM
DF = 0.15 * FCEN

N_GLASS = 1.5
N_ISL = 1.75
N_GAP = 1.0
N_AIR = 1.0


def rt(fill, Lambda=0.50, d=0.20, res=120, dpml=0.6, pad=0.8):
    sy = dpml + pad + d + pad + dpml
    cell = mp.Vector3(Lambda, sy)
    pml = [mp.PML(dpml, direction=mp.Y)]

    y0 = -sy / 2
    y_sub_top = y0 + dpml + pad          # bottom of patterned layer
    y_lay0 = y_sub_top
    y_lay1 = y_lay0 + d
    y_src = y0 + dpml + 0.25 * pad
    y_refl = y0 + dpml + 0.55 * pad
    y_tran = y_lay1 + 0.5 * pad

    # source: plane wave (line spanning the period) in glass, going +y
    src = [mp.Source(mp.GaussianSource(FCEN, fwidth=DF), component=mp.Ez,
                     center=mp.Vector3(0, y_src), size=mp.Vector3(Lambda, 0))]

    def make_geometry():
        top = -y0                                 # top of cell (y coordinate)
        g = [
            # glass substrate up to the patterned layer
            mp.Block(center=mp.Vector3(0, (y0 + y_lay0) / 2),
                     size=mp.Vector3(mp.inf, y_lay0 - y0, mp.inf),
                     material=mp.Medium(index=N_GLASS)),
            # air superstrate above the layer
            mp.Block(center=mp.Vector3(0, (y_lay1 + top) / 2),
                     size=mp.Vector3(mp.inf, top - y_lay1, mp.inf),
                     material=mp.Medium(index=N_AIR)),
            # patterned layer: gap fill, then island block over the fraction
            mp.Block(center=mp.Vector3(0, (y_lay0 + y_lay1) / 2),
                     size=mp.Vector3(mp.inf, d, mp.inf),
                     material=mp.Medium(index=N_GAP)),
            mp.Block(center=mp.Vector3(0, (y_lay0 + y_lay1) / 2),
                     size=mp.Vector3(fill * Lambda, d, mp.inf),
                     material=mp.Medium(index=N_ISL)),
        ]
        return g

    # ---- normalization run: uniform glass everywhere (pure incident) ----
    sim = mp.Simulation(cell_size=cell, resolution=res, boundary_layers=pml,
                        sources=src, k_point=mp.Vector3(),
                        default_material=mp.Medium(index=N_GLASS))
    refl = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
        center=mp.Vector3(0, y_refl), size=mp.Vector3(Lambda, 0)))
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        5, mp.Ez, mp.Vector3(0, y_refl), 1e-7))
    incident = mp.get_fluxes(refl)[0]
    refl_data = sim.get_flux_data(refl)

    # ---- structure run --------------------------------------------------
    sim.reset_meep()
    sim = mp.Simulation(cell_size=cell, resolution=res, boundary_layers=pml,
                        sources=src, k_point=mp.Vector3(),
                        geometry=make_geometry())
    refl = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
        center=mp.Vector3(0, y_refl), size=mp.Vector3(Lambda, 0)))
    tran = sim.add_flux(FCEN, 0, 1, mp.FluxRegion(
        center=mp.Vector3(0, y_tran), size=mp.Vector3(Lambda, 0)))
    sim.load_minus_flux_data(refl, refl_data)   # subtract incident
    sim.run(until_after_sources=mp.stop_when_fields_decayed(
        5, mp.Ez, mp.Vector3(0, y_tran), 1e-7))

    refl_flux = mp.get_fluxes(refl)[0]
    tran_flux = mp.get_fluxes(tran)[0]
    R = -refl_flux / incident
    T = tran_flux / incident
    return R, T


def tmm_planar(d, n2=N_ISL):
    """Analytic normal-incidence R/T for glass | film(n2, d) | air."""
    k0 = 2 * np.pi / LAM
    n1, n3 = N_GLASS, N_AIR

    def r(a, b):
        return (a - b) / (a + b)

    r12, r23 = r(n1, n2), r(n2, n3)
    ph = np.exp(2j * n2 * k0 * d)
    rt_ = (r12 + r23 * ph) / (1 + r12 * r23 * ph)
    R = abs(rt_) ** 2
    return R, 1 - R                       # lossless: T = 1 - R


def tmm_ema(d, fill):
    """'Best-case' planar TMM: treat the island layer as an effective medium
    (volume-averaged permittivity). This is what a TMM user would do to fake a
    patterned layer; it ignores diffraction and so must fail when the period is
    comparable to the wavelength."""
    eps_eff = fill * N_ISL ** 2 + (1 - fill) * N_GAP ** 2
    return tmm_planar(d, n2=np.sqrt(eps_eff))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", type=int, default=120)
    ap.add_argument("--Lambda", type=float, default=0.50)
    ap.add_argument("--d", type=float, default=0.20)
    args = ap.parse_args()

    print("lambda=550nm, normal incidence, glass(1.5)|layer(1.75/air, d=%g um)|air"
          % args.d)
    print("period Lambda = %g um  (lambda/n_glass = %.3f um -> higher orders in "
          "glass if Lambda > that)\n" % (args.Lambda, LAM / N_GLASS))

    Rt, Tt = tmm_planar(args.d)
    print("[TMM planar, continuous slab]      R=%.4f  T=%.4f" % (Rt, Tt))

    Rf, Tf = rt(fill=1.0, Lambda=args.Lambda, d=args.d, res=args.res)
    print("[FDTD fill=1.0  (validation)]       R=%.4f  T=%.4f  -> approaches TMM "
          "(1st-order in res)\n" % (Rf, Tf))

    print("fill |  FDTD R  FDTD T | EMA-TMM R  EMA-TMM T |  |dR| vs EMA")
    for fill in (0.7, 0.5, 0.3):
        Ri, Ti = rt(fill=fill, Lambda=args.Lambda, d=args.d, res=args.res)
        Re, Te = tmm_ema(args.d, fill)
        print(" %.1f | %.4f  %.4f | %.4f    %.4f   |  %.3f" %
              (fill, Ri, Ti, Re, Te, abs(Ri - Re)))


if __name__ == "__main__":
    main()
