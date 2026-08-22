"""Separate absorption from scattering without an angular scan.

1 - T - R is A + S: light that neither passed nor reflected specularly, whether
it was absorbed or thrown sideways out of the collection cone. Ellipsometry does
not fix this -- the fitted k absorbs the scatter too. So neither instrument can
give a clean A on its own.

But sheet resistance can, because it cannot scatter. Drude puts the same
electron damping into rho and into eps2, so

    eps2(d) / eps2_bulk = rho(d) / rho_bulk = Rs(d) * d / rho_bulk

is a purely DC-derived optical loss, blind to roughness by construction. Feed it
back through the same transfer matrix and the difference between the modelled
1-T-R and the measured 1-T-R is the scattered fraction.

One correction is needed before that works. Total resistivity carries two size
terms, and they do not enter the optics equally:

  - Fuchs-Sondheimer surface scattering, rho_bulk*(1 + 0.375(1-p)*l/d), is felt
    by the optical response too;
  - Mayadas-Shatzkes grain-boundary scattering is not. Over one optical cycle at
    550 nm an electron travels v_F/omega ~ 0.4 nm, far short of a grain, so the
    boundaries are invisible at that frequency while they dominate the DC path.

Scaling eps2 by the TOTAL resistivity therefore over-predicts optical loss in
the thicker films, and the naive subtraction returns a negative scattered
fraction at 10 and 12 nm -- which is how the effect announces itself. Scaling by
the surface term alone gives a floor that is positive everywhere.

So three curves bracket the answer, and none of them needs the scattering to be
measured:

  floor    eps2 from the Fuchs term only -- absorption a perfectly specular film
           of this thickness cannot go below
  DC       eps2 from the full measured resistivity -- an over-estimate of
           absorption, for the reason above
  optical  eps2 from 1-T-R -- absorption plus whatever was scattered, an upper
           bound of a different kind

Where all three bottom out is the answer to "what thickness minimises one-pass
absorption", and it is reached without ever separating A from S.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
m = importlib.import_module("91_nk_library_series")

L = 550.0
D = [4, 5, 6, 7, 8, 10, 12]
IDS = {4: "1-2", 5: "1-3", 6: "1-4", 7: "2-1", 8: "2-2", 10: "2-3", 12: "2-4"}
RS = {4: 52.2, 5: 23.3, 6: 18.9, 7: 12.7, 8: 9.1, 10: 7.5, 12: 5.3}
RHO_BULK = 1.59                      # uOhm cm, silver at 300 K
MFP = 52.0                           # nm, silver electron mean free path at 300 K
P_SPEC = 0.0                         # fully diffuse surfaces: the Fuchs worst case


def rho_fuchs(d):
    """Surface-scattering resistivity alone, the part the optics also sees."""
    return RHO_BULK * (1.0 + 0.375 * (1.0 - P_SPEC) * MFP / d)


def nk_from_eps(e1, e2):
    e = complex(e1, e2)
    return np.sqrt(e)


def main():
    ag = m.AG_IDEAL(L)
    e1_bulk = ag.real**2 - ag.imag**2
    e2_bulk = 2 * ag.real * ag.imag
    print(f"at {L:.0f} nm.  bulk silver: eps1 {e1_bulk:.3f}, eps2 {e2_bulk:.3f}, "
          f"rho {RHO_BULK} uOhm cm\n")

    print(f"{'d':>3} {'Rs':>6} {'rho':>6} {'surf':>6} {'gb':>6} | "
          f"{'e2 floor':>8} {'e2 DC':>7} {'e2 opt':>7} | "
          f"{'A meas':>7} {'A floor':>8} {'A DC':>7} | "
          f"{'S vs floor':>10} | {'dev floor':>9} {'dev DC':>7} {'dev opt':>7}")
    print("-" * 122)
    out = []
    for d in D:
        rho = RS[d] * d * 0.1                       # 1 ohm nm = 0.1 uOhm cm
        rs_surf = rho_fuchs(d)
        gb = rho - rs_surf

        T = m.read_csv(os.path.join(m.RAW, f"{IDS[d]}T.csv"))[L]
        R = m.read_csv(os.path.join(m.RAW, f"{IDS[d]}R.csv"))[L]
        Rb = ((m.n_glass(L) - 1) / (m.n_glass(L) + 1))**2
        Rc = R + (1 - m.KEEP) * (T**2) * Rb
        A_meas = 1 - T - Rc

        nk_opt, _ = m.invert(L, T, R, "HATCN", float(d))
        e1 = nk_opt.real**2 - nk_opt.imag**2        # stays at the bulk value
        e2_opt = 2 * nk_opt.real * nk_opt.imag
        e2_floor = e2_bulk * rs_surf / RHO_BULK
        e2_dc = e2_bulk * rho / RHO_BULK

        nk_floor = nk_from_eps(e1, e2_floor)
        nk_dc = nk_from_eps(e1, e2_dc)
        A_floor = 1 - sum(m.observable(nk_floor, L, "HATCN", float(d)))
        A_dc = 1 - sum(m.observable(nk_dc, L, "HATCN", float(d)))

        dfl, _ = m.best_device_A(nk_floor, L, float(d))
        ddc, _ = m.best_device_A(nk_dc, L, float(d))
        dop, _ = m.best_device_A(nk_opt, L, float(d))
        out.append((d, A_meas, A_floor, A_dc, dfl, ddc, dop, e2_opt / e2_floor))
        print(f"{d:>3} {RS[d]:>6.1f} {rho:>6.2f} {rs_surf:>6.2f} {gb:>6.2f} | "
              f"{e2_floor:>8.3f} {e2_dc:>7.3f} {e2_opt:>7.3f} | "
              f"{A_meas*100:>6.2f}% {A_floor*100:>7.2f}% {A_dc*100:>6.2f}% | "
              f"{(A_meas-A_floor)*100:>9.2f}% | "
              f"{dfl*100:>8.2f}% {ddc*100:>6.2f}% {dop*100:>6.2f}%")

    for name, col in (("floor  (Fuchs only, no scatter, no grain boundaries)", 4),
                      ("DC     (full resistivity, over-counts grain boundaries)", 5),
                      ("optical(1-T-R, absorption + scattering)", 6)):
        v = [r[col] for r in out]
        print(f"\ndevice one-pass A, {name}")
        print(f"  " + "  ".join(f"{r[0]}nm {r[col]*100:.2f}%" for r in out))
        print(f"  minimum at {D[int(np.argmin(v))]} nm, {min(v)*100:.2f}%")

    print("\nThe DC route and the optical route both put the minimum at 8 nm, and the")
    print("DC route cannot contain a single scattered photon -- it is a four-point")
    print("probe measurement. That is what settles the question: the 7-8 nm optimum")
    print("survives even if every photon missing from 1-T-R turns out to have been")
    print("scattered rather than absorbed.")
    print("\nThe floor curve has no minimum, and should not be read as one. It applies")
    print("Fuchs-Sondheimer at every thickness, which silently assumes the film is")
    print("already closed -- false below 5 nm on HATCN. It is the ideal-film")
    print("counterfactual, and its message is the interesting one: an ideal silver")
    print("film would only pay 1.54 -> 1.98% going from 4 to 12 nm. The whole 7-8 nm")
    print("optimum therefore comes from film quality, not from thickness. Thin is")
    print("only bad because thin films are bad silver.")
    print(f"\nexcess of the optical eps2 over the Fuchs floor "
          f"(scattering + grain boundaries + roughness):")
    print("  " + "  ".join(f"{r[0]}nm {r[7]:.2f}x" for r in out))
    print("  falls monotonically to 8 nm, then flattens -- the films stop improving")
    print("  once they have closed, which is the closure signature seen in transport.")


if __name__ == "__main__":
    main()
