"""How large is interfacial CT absorption, and how would you measure it?

WHY THIS EXISTS. scripts/47 computes a bulk k for neutral HATCN and shows the
visible tail is negligible. It also says, twice, that this does NOT bound the
absorption from charge-transfer states at the HATCN/donor and HATCN/Ag interfaces,
because that absorption scales with interface AREA, not with layer thickness, and
therefore does not shrink as the seed is made thinner. That is the loss channel
that actually applies to a HIL, and leaving it as a caveat is not good enough:
without a magnitude the manuscript cannot say whether it matters.

THE PHYSICS. Where HATCN (acceptor) meets a donor HTL or an Ag surface, the ground
state acquires partial CT character and a new optical transition appears BELOW the
gap of either component -- typically in the visible or NIR. This is not a defect:
it is the same charge transfer that makes HATCN function as a hole-injection layer.
The transition is localised to roughly one molecular layer at the junction.

WHAT IS ESTIMATED HERE. A monolayer of CT chromophores has a 2D density set by the
molecular footprint, so its absorptance is

    A = sigma(nu) * N_2D

with sigma from the oscillator strength. Sweeping f over the plausible range gives
the size of the effect per interface, which is what sets the sensitivity an
experiment has to reach.

WHAT IS NOT ESTIMATED. The oscillator strength of the CT transition itself. It
depends on donor-acceptor orbital overlap and geometry, is not transferable between
systems, and computing it needs an explicit interface calculation (a donor/HATCN
or Ag-slab/HATCN dimer at TD-DFT, which is a separate job). So this is a
SENSITIVITY MAP over f, not a prediction. It answers "how well must I measure",
not "what will I see".
"""
import numpy as np

E_CHG, EPS0, M_E, C0 = 1.602176634e-19, 8.8541878128e-12, 9.1093837015e-31, 2.99792458e8
H_EVS = 4.135667696e-15

# HATCN lies flat on both metals and organics; its van der Waals disc is ~1.1-1.3 nm
# across, giving a footprint near 1.2 nm^2. Face-on packing is the dense case and
# therefore the pessimistic one for absorption.
FOOTPRINT_NM2 = 1.2
FWHM_EV = 0.5          # CT bands are broad and structureless
F_GRID = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]


def monolayer_absorptance(f, fwhm_ev=FWHM_EV, footprint_nm2=FOOTPRINT_NM2):
    """Peak absorptance of one monolayer of CT chromophores, single pass."""
    n2d = 1.0 / (footprint_nm2 * 1e-18)                 # molecules / m^2
    integ_nu = f * E_CHG ** 2 / (4 * EPS0 * M_E * C0)   # m^2 * Hz
    integ_ev = integ_nu * H_EVS                          # m^2 * eV
    sigma_pk = integ_ev / (fwhm_ev / (2 * np.sqrt(2 * np.log(2))) * np.sqrt(2 * np.pi))
    return sigma_pk * n2d


def main():
    print("Peak absorptance of ONE CT monolayer (single pass, per interface)\n")
    print(f"  footprint {FOOTPRINT_NM2} nm^2/molecule -> "
          f"{1.0/(FOOTPRINT_NM2*1e-18):.2e} m^-2,  FWHM {FWHM_EV} eV\n")
    print(f"{'f (CT)':>9}{'A per interface':>18}{'x2 interfaces':>16}   note")
    print("-" * 62)
    for f in F_GRID:
        a = monolayer_absorptance(f)
        note = ("below UV-Vis reproducibility" if 100 * a < 0.1 else
                "measurable" if 100 * a < 1.0 else "large -- would show in T")
        print(f"{f:>9.3f}{100*a:>17.3f}%{200*a:>15.3f}%   {note}")

    print("\n" + "=" * 62)
    print("HOW TO SEPARATE INTERFACE FROM BULK")
    print("=" * 62)
    print("""
The two scale differently, and that is the whole lever:

    A_total(d, N) = alpha_bulk * d  +  N * A_interface

  (1) THICKNESS SERIES, fixed interfaces. Glass/HATCN(d)/Ag with
      d = 2, 5, 10, 20, 30 nm. Plot A vs d: the SLOPE is the bulk k, the
      INTERCEPT at d -> 0 is the interfacial term. This is the cleanest
      measurement and needs only equipment already in use for the film series.

  (2) MULTILAYER AMPLIFICATION, for the organic/organic junction. Compare
      [HATCN(5)/NPB(5)] x N against a single HATCN(5N)/NPB(5N) bilayer. Same
      total material, N times the interface area, so the difference is
      (N-1) interfaces. N = 10 turns a 0.1 % effect into ~1 %, which moves it
      from undetectable to comfortable. This is the way to measure it if the
      single-interface signal sits under the noise floor.

  (3) SPECTRAL WINDOW. The CT band lies BELOW the gap of both components, so
      look where neither absorbs. For HATCN (edge ~380-400 nm) against a typical
      HTL, that is roughly 450-900 nm -- a clean window, which is why a weak band
      is still findable. Subtract: A(bilayer) - A(HATCN) - A(donor).

  (4) THE Ag INTERFACE IS HARDER. Silver reflects and absorbs strongly, so a
      monolayer's contribution is a small residual on a large background, and
      ellipsometric fitting of an explicit interface layer is badly correlated
      with the Ag optical constants. The right tool is DIFFERENTIAL REFLECTANCE
      (dR/R measured in situ as HATCN is deposited on the Ag film): interfacial
      CT saturates at ~1 monolayer while any bulk contribution keeps growing
      linearly, so the COVERAGE DEPENDENCE separates them without needing an
      absolute reference. If in-situ optics are not available, fall back to (1).

  (5) MECHANISM CHECK, not a magnitude. UPS on the interface shows whether CT
      induced gap states exist and where. Pair it with the optical result --
      neither alone distinguishes a CT band from a scattering loss.
""")
    print("=" * 62)
    print("FOR THE MANUSCRIPT")
    print("=" * 62)
    print("""
The honest operational statement is the TMM residual. Predict the stack
transmittance from bulk n,k alone; measure it; report the difference. That
residual bounds interfacial CT plus scattering plus model error TOGETHER, and
does not require attributing it. If the residual is within reproducibility, the
claim is "any interfacial contribution is below X %p", which is exactly as strong
as the data allow and is not vulnerable to a reviewer who knows that CT states at
an acceptor/metal junction are unavoidable.

Do NOT write that HATCN contributes no absorption. It contributes none from its
bulk (scripts/47) and an unmeasured amount from its interfaces, and those are
different statements.
""")


if __name__ == "__main__":
    main()
