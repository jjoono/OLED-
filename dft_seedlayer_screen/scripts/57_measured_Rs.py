"""Measured sheet resistance on HATCN(4 nm)/Ag: the first real test of the model.

THE DATA (measured, HATCN 4 nm seed):

    Ag  3 nm  ->  ~100 Ohm/sq
    Ag  5 nm  ->  ~25  Ohm/sq
    Ag 25 nm  ->  ~1.5 Ohm/sq

WHAT IT SETTLES IMMEDIATELY. All three conduct. A film that had not percolated
would read open circuit or megaohms, so Ag is CONTINUOUS at 3 nm on a 4 nm HATCN
seed. That is the single most important number here and it is better than the SEM
alone could establish -- the granular texture at 3 nm is a just-percolated film,
not a disconnected one.

WHAT IT CORRECTS. scripts/37 sets D_nm = 100 for the Mayadas-Shatzkes grain size,
a value taken from SEM of thicker films after the textbook D = d appeared to
overestimate Rs at 20-25 nm. Against this measurement D = 100 nm is wrong and
wrong in a thickness-dependent way -- it underestimates Rs by 4.7x at 3 nm, 2.6x
at 5 nm, 1.35x at 25 nm. Restoring D = d fits all three within ~30 % with no
fitted parameter. The earlier correction is retracted for this system.

WHAT REMAINS OPEN. Inverting the measurement for D gives 2.0 / 6.0 / 33.7 nm at
3 / 5 / 25 nm. The 25 nm value disagrees with the ~100 nm features seen by SEM by
a factor of three. Both cannot be the grain size. The likely resolution is that
SEM shows surface mounds or agglomerates that each contain several grains, since
grain boundaries are not generally visible in secondary-electron imaging -- but
that is a hypothesis, and EBSD or a TEM plan view would settle it.
"""
import importlib.util, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("m", os.path.join(HERE, "37_ag_sheet_resistance.py"))
M = importlib.util.module_from_spec(_s); _s.loader.exec_module(M)

RHO0 = 1.587e-8            # Ohm.m, bulk Ag
MEAS = {3.0: 100.0, 5.0: 25.0, 25.0: 1.5}
SEED = "HATCN 4 nm"

# Ag optical constants, Johnson & Christy, as used throughout the project
JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
WL = np.arange(400.0, 701.0, 2.0)


def tmm_T(d_ag, n_sup=1.0, n_sub=1.52):
    """Normal-incidence T of a bare Ag film. Macleod convention: conjugate n+ik."""
    n_ag = np.conj(np.interp(WL, JC_WL, AG_N) + 1j * np.interp(WL, JC_WL, AG_K))
    T = np.zeros_like(WL)
    for i, lam in enumerate(WL):
        d = 2 * np.pi / lam * n_ag[i] * d_ag
        c, s = np.cos(d), np.sin(d)
        Mx = np.array([[c, 1j * s / n_ag[i]], [1j * n_ag[i] * s, c]])
        B, C = Mx @ np.array([1.0, n_sub], dtype=complex)
        T[i] = 4 * n_sup * n_sub / abs(n_sup * B + C) ** 2
    v = 1.019 * np.exp(-285.4 * ((WL / 1000.0) - 0.559) ** 2)
    return float((T * v).sum() / v.sum())


def main():
    print(f"MEASURED, seed = {SEED}\n")
    print(f"{'d (nm)':>7}{'Rs meas':>10}{'rho':>12}{'rho/rho0':>10}"
          f"{'Rs bulk':>10}")
    for d, rs in MEAS.items():
        rho = rs * d * 1e-9
        print(f"{d:>7.0f}{rs:>10.1f}{rho*1e8:>10.2f} uOcm{rho/RHO0:>10.1f}"
              f"{RHO0/(d*1e-9):>10.2f}")

    print("\n" + "=" * 62)
    print("MODEL COMPARISON (Fuchs-Sondheimer + Mayadas-Shatzkes)")
    print("=" * 62)
    print(f"{'d':>4}{'meas':>8}{'D=100nm':>10}{'ratio':>7}{'D=d':>9}{'ratio':>7}")
    for d, rs in MEAS.items():
        r100 = M.combined(d, p=0.5, Rgb=0.3, D_nm=100.0) / (d * 1e-9)
        rd = M.combined(d, p=0.5, Rgb=0.3, D_nm=d) / (d * 1e-9)
        print(f"{d:>4.0f}{rs:>8.1f}{r100:>10.1f}{rs/r100:>7.2f}"
              f"{rd:>9.1f}{rs/rd:>7.2f}")
    print("""
  D = d fits all three within ~30 % with nothing fitted; D = 100 nm fails
  systematically and worst where the film is thinnest. scripts/37's default
  should go back to D = d for this system.

  The residual at 3 nm -- measurement 100 vs model 73 -- runs the right way for
  a just-percolated film: current follows a tortuous path through a barely
  connected network, and FS+MS assumes a uniform continuous sheet, so it cannot
  see that. The agreement at 5 and 25 nm, where the film is properly continuous,
  is the part that tests the model.""")

    print("\n" + "=" * 62)
    print("FIGURE OF MERIT")
    print("=" * 62)
    print(f"{'d':>4}{'Rs':>8}{'T (photopic)':>15}{'sigma_dc/sigma_opt':>21}")
    for d, rs in MEAS.items():
        T = tmm_T(d)
        # Haacke-style ratio used for transparent conductors:
        #   sigma_dc/sigma_opt = 188.5 / (Rs * (T^-0.5 - 1))
        fom = 188.5 / (rs * (T ** -0.5 - 1)) if T > 0 else float("nan")
        print(f"{d:>4.0f}{rs:>8.1f}{100*T:>14.1f}%{fom:>21.1f}")
    print("""
  T here is the BARE Ag film on glass with no capping layer, so it is the
  pessimistic case -- scripts/32 showed a 40 nm cap lifts relative transmittance
  to ~100 %. Treat these as a floor, not as the device numbers.

  For orientation: ITO runs 10-100 Ohm/sq at 80-90 % transmittance. The 5 nm
  point is in that range on both axes, which is what makes it the interesting
  thickness rather than 3 nm.""")

    print("\n" + "=" * 62)
    print("WHAT TO MEASURE NEXT")
    print("=" * 62)
    print("""
  1. The SAME series WITHOUT the HATCN seed. Every claim that the seed helps
     rests on a comparison that does not exist yet. Bare Ag on glass typically
     does not percolate until 10-15 nm, so if 3 nm conducts here and not there,
     that is the result -- and it is a strong one.
  2. Ag 2 nm and 2.5 nm, to find where conduction actually stops. The percolation
     threshold is now known only to be below 3 nm.
  3. Grain size by EBSD or TEM plan view at 25 nm, to resolve 34 nm (transport)
     against ~100 nm (SEM).
  4. T of the same samples, so the figures of merit above use measured rather
     than modelled transmittance.""")


if __name__ == "__main__":
    main()
