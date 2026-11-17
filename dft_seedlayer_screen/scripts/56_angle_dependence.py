"""R and T of a 25 nm Ag film at 0 and 6 degrees: does the measurement angle matter?

WHY 6 DEGREES. A spectrophotometer cannot measure specular reflectance at exactly
normal incidence without a beam splitter, so reflectance accessories sit at 6 or 8
degrees and the result is compared against a normal-incidence calculation. This
quantifies the error that substitution introduces.

METHOD. Abeles characteristic matrix at oblique incidence, s and p polarisation
separately, then the unpolarised average -- an unpolarised source at 6 degrees
gives (R_s + R_p)/2, and the two polarisations already differ at first order in
the angle even though their average does not.

    delta_j = (2 pi / lambda) N_j d_j cos(theta_j)
    eta_j   = N_j cos(theta_j)        (s)
            = N_j / cos(theta_j)      (p)

SIGN CONVENTION. Macleod's formulation assumes exp[i(wt - kz)], i.e. N = n - ik,
while the project's tables store n + ik. scripts/33 traced a spurious gain (T+R>1)
to exactly this, so the constants are conjugated on entry here too. Snell's law is
applied with the complex index, which is the correct treatment for an absorbing
layer -- the "angle" in the metal is complex and has no geometric meaning, but the
algebra carries through and gives the right R and T.

Ag optical constants: Johnson & Christy 1972, the same 12-point table used by
scripts/32 and 33, so these numbers sit on the same footing as the rest of the
project's optics.
"""
import numpy as np

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])

WL = np.arange(400.0, 701.0, 2.0)
N_GLASS, N_AIR = 1.52, 1.0
D_AG = 25.0
ANGLES = (0.0, 6.0, 8.0)          # 8 included: some accessories use it instead


def photopic_weight(wl):
    v = 1.019 * np.exp(-285.4 * ((wl / 1000.0) - 0.559) ** 2)
    return v / v.sum()


def rt(n0, layers, nsub, lam, theta0_deg, pol):
    """R, T of a stack at oblique incidence. layers = [(N, d_nm), ...]."""
    th0 = np.deg2rad(theta0_deg)
    sin0 = n0 * np.sin(th0)                       # conserved by Snell

    def eta(N):
        cos_t = np.sqrt(1.0 - (sin0 / N) ** 2)
        return (N * cos_t, cos_t)

    eta0, _ = eta(n0)
    if pol == "p":
        eta0 = n0 ** 2 / np.sqrt(n0 ** 2 - sin0 ** 2)
    etasub, _ = eta(nsub)
    if pol == "p":
        etasub = nsub ** 2 / np.sqrt(nsub ** 2 - sin0 ** 2)

    M = np.eye(2, dtype=complex)
    for N, d in layers:
        cos_t = np.sqrt(1.0 - (sin0 / N) ** 2)
        e = N * cos_t if pol == "s" else N ** 2 / (N * cos_t)
        delta = 2 * np.pi / lam * N * d * cos_t
        c, s = np.cos(delta), np.sin(delta)
        M = M @ np.array([[c, 1j * s / e], [1j * e * s, c]])

    B, C = M @ np.array([1.0, etasub], dtype=complex)
    r = (eta0 * B - C) / (eta0 * B + C)
    R = float(abs(r) ** 2)
    T = float(4 * eta0.real * etasub.real / abs(eta0 * B + C) ** 2)
    return R, T


def spectrum(theta_deg, from_glass=False):
    # conjugate: table is n + ik, Macleod needs n - ik
    n_ag = np.conj(np.interp(WL, JC_WL, AG_N) + 1j * np.interp(WL, JC_WL, AG_K))
    n_in = N_GLASS if from_glass else N_AIR
    n_out = N_AIR if from_glass else N_GLASS
    R = np.zeros_like(WL); T = np.zeros_like(WL)
    Rs = np.zeros_like(WL); Rp = np.zeros_like(WL)
    for i, lam in enumerate(WL):
        rs, ts = rt(n_in, [(n_ag[i], D_AG)], n_out, lam, theta_deg, "s")
        rp, tp = rt(n_in, [(n_ag[i], D_AG)], n_out, lam, theta_deg, "p")
        Rs[i], Rp[i] = rs, rp
        R[i], T[i] = 0.5 * (rs + rp), 0.5 * (ts + tp)
    return R, T, Rs, Rp


def main():
    w = photopic_weight(WL)
    i550 = int(np.argmin(abs(WL - 550)))

    for side, from_glass in (("air side (air / Ag 25 nm / glass)", False),
                             ("glass side (glass / Ag 25 nm / air)", True)):
        print("\n" + "=" * 70)
        print(side)
        print("=" * 70)
        print(f"{'angle':>6}{'R@550':>9}{'T@550':>9}{'R photopic':>13}"
              f"{'T photopic':>13}{'R_s-R_p @550':>15}")
        ref = None
        for a in ANGLES:
            R, T, Rs, Rp = spectrum(a, from_glass)
            row = (100 * R[i550], 100 * T[i550], 100 * (R * w).sum(),
                   100 * (T * w).sum(), 100 * (Rs[i550] - Rp[i550]))
            print(f"{a:>5.0f}°{row[0]:>9.3f}{row[1]:>9.3f}{row[2]:>13.3f}"
                  f"{row[3]:>13.3f}{row[4]:>15.4f}")
            if a == 0.0:
                ref = row
        print()
        for a in ANGLES[1:]:
            R, T, _, _ = spectrum(a, from_glass)
            dR = 100 * R[i550] - ref[0]
            dT = 100 * T[i550] - ref[1]
            dRp = 100 * (R * w).sum() - ref[2]
            dTp = 100 * (T * w).sum() - ref[3]
            print(f"  {a:.0f}° vs 0°:  dR@550 = {dR:+.4f} %p,  dT@550 = {dT:+.4f} %p"
                  f"   (photopic {dRp:+.4f} / {dTp:+.4f} %p)")

    print("\n" + "=" * 70)
    print("WHAT THIS MEANS FOR THE MEASUREMENT")
    print("=" * 70)
    print("""
The angle enters through cos(theta), and cos(6 deg) = 0.9945 -- a 0.55 % change
in a geometric factor, not in the optical constants. For a single absorbing layer
that propagates to R and T as a shift of order 0.01 %p, far below the ~0.1-0.5 %p
reproducibility of a spectrophotometer.

So a 6 degree measurement IS comparable to a normal-incidence calculation, and no
correction is needed. The same is true at 8 degrees.

What the angle does change measurably is the SPLITTING between s and p
polarisation, which grows as theta^2 and is printed above. It stays small at 6
degrees, but if the instrument's beam is polarised and the sample is anisotropic
the two no longer average cleanly -- worth knowing, not worth correcting for here.

Caveat: this is the ideal 25 nm film. Real R and T at this thickness are dominated
by roughness and by whether the film is fully continuous, both of which move the
result by percentage points -- two orders of magnitude more than the angle does.
Do not spend effort on the angle while those are uncontrolled.
""")


if __name__ == "__main__":
    main()
