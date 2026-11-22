# Reciprocity separates reflective and transmissive angular mixing in OLED light extraction

**Authors** — [TBD]

---

## Abstract

Angular mixing is widely assumed to help extract the light trapped in the substrate of organic light-emitting diodes, as in the fully randomized model of statistical ray optics. Using optical reciprocity, we show that the benefit depends on where the mixing takes place. An escape interface that mixes light in reflection returns it preferentially to angles from which it cannot escape, which lowers the escape probability by an amount set by its angular variance. Mixing in transmission through a weakly absorbing volume carries no such penalty and recovers the fully randomized model. For 37 microlens, scattering, and textured structures, complete mixing by the structure adds at most about 2 percentage points to the spectrally weighted extraction, whereas a forward-scattering volume recovers 92–99% of the fully randomized gain over a planar interface. A one-parameter reduced model reproduces 342 angle-resolved calculations with a mean error of 0.30 percentage points, and microlens devices with two emitters confirm the predicted angular redistribution.

---

## 1. Introduction

A planar organic light-emitting diode (OLED) emits only a fraction of the light it generates. The rest is lost to waveguided and surface-plasmon modes, to absorption in the electrodes, and to total internal reflection at the substrate–air interface [1,2]. The last of these losses can be reduced without touching the electrical stack by structuring the substrate with microlens arrays (MLAs) [3,4], scattering layers [5,6], or textured films [7]. The dependence of extraction on lens profile, aspect ratio, and fill factor is well documented [3,8], and deep-learning methods are now applied to the design of both the thin-film stack [9] and the outcoupling microstructure [10]. Even so, the external quantum efficiencies (EQEs) reported with such structures have mostly stayed between about 50% and 60% [5,11–13], a ceiling that a recent analysis attributes to the loss incurred each time unextracted light returns to the device [14].

Two kinds of model describe this multi-pass process. Statistical ray optics assumes that the light inside the substrate is fully randomized at every pass [15–17], so that the substrate-to-air efficiency, $\eta = p/[p+(1-p)A']$, depends only on the escape probability $p$ of light with uniform radiance and on the round-trip loss $A'$ [14]. Angle-resolved models instead carry the angular power distribution from one pass to the next with redistribution matrices obtained from ray tracing or radiative transfer [2,14,18–22]. They are accurate, but they return a number for each structure and no rule for comparing structures. How far a real texture is from the randomized limit has been addressed mainly case by case, through the absorption enhancement of textured films [23], pass-resolved ray tracing of textured silicon [24], and the conditions under which a rough interface scatters like a Lambertian surface [25].

The two descriptions disagree most where it matters for design. A hemispherical MLA returns most unextracted light near the angle at which it arrived [14]. The light that survives several passes therefore sits at angles the lenses transmit poorly, and the fully randomized model overestimates the extraction by an amount that differs from structure to structure. It is tempting to conclude that outcoupling structures should be designed to mix angles faster. Whether faster mixing pays off, and by how much, has not been established.

Here we show that the answer depends on where the light is mixed, and that optical reciprocity fixes this dependence. Reciprocity makes the return matrix of a passive structure self-adjoint with respect to étendue. An escape interface that redistributes light in reflection must then return it to each angle in proportion to the fraction that this angle reflects, and light that is completely mixed in this way escapes with probability $p - \sigma_p^2/(1-p)$, where $\sigma_p^2$ is the angular variance of the escape probability. Mixing in transmission through a lossless volume has no such penalty, and mixing at the device is penalized only through the covariance between its reflectance and the escape probability. The fully randomized model is thus a model of transmissive mixing. It does not describe the reflective mixing that refractive outcoupling structures provide.

Section 2 derives these results together with a reduced model that has a single mixing parameter. The second eigenvalue of the return matrix sets this parameter and defines an angular mixing number, the per-pass counterpart of the transport factor $1/(1-g)$ of radiative transfer [26] and of the relaxation time of a reversible Markov chain [27]. The leading eigenvalue gives the escape probability of light that has circulated many times, which is the quasi-stationary distribution of an absorbing chain [28] and corresponds to the escape rate of an open dynamical system [29–31]. Spectral descriptions of this kind are established for light-transport operators [32–34]; their consequences for light extraction have not been worked out. Section 4 evaluates 38 structures in three device stacks, tests the reduced model against full transport, direct ray tracing, and structures chosen after the model had been fixed, and compares the finite-geometry form of the model with fabricated MLA devices of two emitters. Section 5 discusses the étendue limit and the implications for design.

---

## 2. Theory

### 2.1 Multi-pass extraction as a Markov chain on angle

We divide the polar angle inside the substrate into 1° bins $k$. Let $v_k$ be the azimuth-integrated power traveling toward the structure in bin $k$, and $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$ the étendue of the bin, so that uniform radiance corresponds to $v_k \propto w_k$. At each encounter the structure transmits a fraction $p_k$ of the light in bin $k$ to air and returns a fraction $B_R(i\to r)$ from bin $i$ to bin $r$ [Fig. 1(a)]. Writing the returned part as

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}) \tag{1}
$$

separates how much light is returned from where it goes. The matrix $M$ is column stochastic and is the transition matrix of a Markov chain on angle. The étendue-weighted mean $p = \sum_k w_k p_k/\sum_k w_k$ is the escape probability of light with uniform radiance. The device reflects the returned light with its angle-resolved reflectance $R(\theta)$ and, unless stated otherwise, preserves the angle. Summing the geometric series over passes gives the substrate-to-air extraction efficiency

$$
\eta_{\mathrm{sta}} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} , \tag{2}
$$

where $\mathbf{v}^{(0)}$ is the power launched into the substrate [14]. All losses enter through $R(\theta)$. The scalar $A' = 1-\sum_k w_k R_k/\sum_k w_k$ is used only for comparison with the fully randomized model.

This description is scalar and averages over azimuth and over the position at which a ray meets the structure. Both simplifications are tested below. Direct tracing in $(\theta,\phi)$ agrees with Eq. (2) to within 0.6 percentage points (pp) (Sec. 4.3), and separate s and p channels change the conclusions by at most 0.4 pp (Table S7).

### 2.2 Reciprocity and the angular mixing number

A passive, reciprocal structure obeys detailed balance in étendue [35],

$$
w_i B_R(i\to j) = w_j B_R(j\to i). \tag{3}
$$

Equation (3) is the discrete form of the self-adjointness of reciprocal scattering operators in the étendue inner product [32]. Redistribution matrices of this kind are used for textured optical sheets [19]. It follows that $M$ is similar to a symmetric matrix, so its eigenvalues are real [27,33]. On the full angular domain its stationary distribution is proportional to $w(1-p)$. For the truncated chain used below it is proportional to $ws$, where $s_i$ is the fraction returned into the retained bins (Note S1). The eigenvalue of second-largest modulus, $\lambda_2$, sets the rate at which an arbitrary angular distribution relaxes toward the stationary one. We express it as a number of encounters,

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|} , \tag{4}
$$

and call $N_{\mathrm{mix}}$ the angular mixing number. A planar interface has $N_{\mathrm{mix}} = \infty$, and a Lambertian scatterer has $N_{\mathrm{mix}} = 1$. Figures 1(b)–1(d) show $M$ for a planar interface, a hemispherical MLA, and a scattering layer. In all structures studied here $\lambda_2$ is real and positive.

Bins close to grazing incidence carry almost no étendue, yet they can keep their angle for many encounters and would otherwise dominate $\lambda_2$. The spectral quantities are therefore evaluated for $\theta < 85°$, which omits 0.76% of the étendue, whereas $\eta_{\mathrm{sta}}$ is always computed over 0–90°. Moving the cutoff between 80° and 90° changes the mean error of the reduced model by less than 0.3 pp (Note S1).

**Fig. 1.** Multi-pass extraction as a Markov chain on angle. (a) A passive structure transmits a fraction $p(\theta)$ of the incident light and returns the rest with $B_R(\theta_i\to\theta_r)$; the device reflects the returned light with $R(\theta)$. (b)–(d) Return matrices $M$ ($\log_{10}$ scale, floor $10^{-3}$) of a planar interface, the hemispherical MLA, and a scattering layer ($g$ = 0.8, $\tau^*$ = 1). The off-diagonal bands in (c) are deterministic paths of rays that are totally reflected at the lens flank, for example from 70° to 34°; they persist in a 6 × 6 supercell and with twice the number of rays (Note S10). (e) Angular distribution of the power incident on the hemispherical MLA at pass $k$, normalized to unit area (Al device, 550 nm). (f) Escape probability at pass $k$ relative to $p$ for five structures. Open symbols on the right mark $p_{\mathrm{mix}}/p$, and the dotted line marks the fully randomized model. *(`fig1_mixing_concept.png`)*

### 2.3 Reciprocal reduced model

The reduced model retains the $p(\theta)$ of a structure and replaces $M$ by a one-parameter family. Of the returned light, a fraction $1-m$ keeps its angle and a fraction $m$ is redistributed over a fixed distribution $q$:

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right]. \tag{5}
$$

Equation (3) is satisfied only for $q \propto w(1-p)$. The light redistributed into an angle is then proportional to the fraction that this angle returns. The same rank-one coupling appears in reciprocal reflectance models that combine a diffuse lobe with an angle-dependent specular term [36,37]. Multi-pass extraction models, in contrast, assume isotropic internal radiance or redistribute the returned light isotropically [15–17,38], and we are not aware of one that imposes the reciprocal form. Isotropic redistribution places part of the returned light inside the escape cone even though the escape cone itself returns almost nothing, which violates Eq. (3) whenever $p$ depends on angle. The normalized kernel of Eq. (5) has the eigenvalues 1 and $1-m$, which identifies $m = 1-\lambda_2$.

### 2.4 Complete mixing at the escape interface

Setting $m = 1$ gives a return matrix of rank one, $B_R(i\to r) = (1-p_i)\,q_r$, and the light it returns escapes at the next encounter with probability

$$
p_{\mathrm{mix}} = \sum_r q_r p_r = p - \frac{\sigma_p^2}{1-p}, \qquad \sigma_p^2 = \frac{\sum_k w_k (p_k - p)^2}{\sum_k w_k}. \tag{6}
$$

The multi-pass series then has the closed form

$$
\eta_{\mathrm{mix}} = \mathbf{p}\cdot\mathbf{v}^{(0)} + \frac{(\mathbf{p}\cdot\mathbf{u})\,[(\mathbf{1}-\mathbf{p})\cdot\mathbf{v}^{(0)}]}{1-(\mathbf{1}-\mathbf{p})\cdot\mathbf{u}}, \qquad u_k = R_k\,q_k . \tag{7}
$$

For an angle-independent reflectance $1-A'$, Eq. (7) becomes $\eta_{\mathrm{mix}} = p_0 + (1-p_0)\,p_{\mathrm{mix}}(1-A')/[p_{\mathrm{mix}} + (1-p_{\mathrm{mix}})A']$, with $p_0$ the first-pass escape probability. It reduces to the fully randomized model only if $p(\theta)$ is constant.

The variance term is large for refractive surfaces, because they transmit inside the critical cone and reflect everything outside it. A planar interface has $p$ = 0.40 but $p_{\mathrm{mix}}$ = 0.06. The fully randomized model assigns the probability $p$ to mixed light and so overestimates what mixing can gain by an amount fixed by $\sigma_p^2$.

Equation (7) is a reference value for rank-one reciprocal redistribution. It is not an upper bound for reciprocal structures with a given $p(\theta)$: Note S1 gives a two-state example in which a kernel of higher rank exceeds $\eta_{\mathrm{mix}}$ by 5 pp. Because $p_{\mathrm{mix}}$ applies at every later pass, the late-pass escape ratio defined in Sec. 2.6 is $\xi = 1 - \sigma_p^2/[p(1-p)]$ for complete mixing and stays below 1 unless $p(\theta)$ is constant.

### 2.5 Mixing in transmission and at the device

The penalty in Eq. (6) arises because the element that mixes the light is the escape interface itself. Reciprocity weights redistribution in reflection by the return probability of the reflecting element, which here is $1-p$. Mixing elsewhere is weighted differently.

Consider first a lossless layer with negligible reflection that mixes the light in transmission, an idealized volume randomizer between the device and the escape interface. Every pass through it distributes the light in proportion to $w$, independently of $p$ and $R$. The mixed light escapes with probability $p$, and complete mixing of this kind reproduces the fully randomized model exactly. In a lossless, index-matched substrate the position of the layer does not matter (Note S11). Real scattering layers also reflect, which lowers their own escape probability, and Sec. 4.1 examines how close they come to the ideal.

The device can also mix the light in reflection, for instance at a diffusely reflecting electrode. The structure keeps its return matrix $B_R$, and the device reflection becomes a kernel $D$ whose column sums equal the reflectance $R_i$ and which must satisfy $w_i D(i\to j) = w_j D(j\to i)$. The rank-one kernel that meets both conditions is

$$
D(i\to j) = \frac{R_i\,R_j\,w_j}{\sum_k w_k R_k} . \tag{8}
$$

Light mixed by the device is distributed in proportion to $wR$ and escapes at the next encounter with probability

$$
p_D = \frac{\sum_k w_k R_k\,p_k}{\sum_k w_k R_k} = p + \frac{\mathrm{Cov}_w(R,p)}{\langle R\rangle_w}, \tag{9}
$$

where the covariance and the mean are étendue weighted. No variance term appears. If $R$ does not depend on angle, $p_D = p$. If the device absorbs most strongly inside the escape cone, as strong microcavities do, the covariance is negative and mixing at the device is penalized as well (Note S11). Table 1 collects the three cases.

**Table 1.** Escape probability of completely mixed light for three placements of the mixing. The last column gives $\eta_{\mathrm{sta}}$ of a planar interface (Al device, 550 nm), which is 0.384 without mixing and 0.798 in the fully randomized model.

| Where the light is mixed | Redistributed in proportion to | Escape probability of the mixed light | Equals the fully randomized model if | $\eta_{\mathrm{sta}}$, planar |
|---|---|---|---|---|
| Escape interface, in reflection | $w(1-p)$ | $p - \sigma_p^2/(1-p)$ | $p(\theta)$ is constant | 0.488 |
| Lossless volume, in transmission | $w$ | $p$ | always | 0.798 |
| Device, in reflection | $wR$ | $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$ | $R$ is uncorrelated with $p$ | 0.770 |

### 2.6 Late-pass escape ratio

The mixing number says how fast the angular distribution relaxes but not what it relaxes to. Light that has circulated many times approaches the leading eigenvector of $B_R$, the quasi-stationary distribution of the absorbing chain [28], and escapes with probability $p_\infty = 1-\rho(B_R)$ per encounter, where $\rho$ denotes the spectral radius. This quantity is the escape rate of an open system [29–31]. The self-adjointness of $B_R$, the Rayleigh quotient evaluated with the uniform-radiance vector, and the column sums of $B_R$ [39] give

$$
\min_\theta p(\theta) \le p_\infty \le \tilde p , \tag{10}
$$

where $\tilde p$ is the probability that light of uniform radiance leaves the truncated chain, either by escape or by return beyond the cutoff ($\tilde p - p < 0.008$). We define the late-pass escape ratio $\xi = p_\infty/\tilde p$. It equals 1 when long-lived light escapes as readily as Lambertian light and tends to 0 when some angles neither escape nor mix. All structures satisfy Eq. (10) numerically (Note S1). $p_\infty$ is a property of the structure alone, for lossless return by the device. In an OLED the surviving light is governed by $\mathrm{diag}\,\mathbf{R}\,B_R$, which also contains the parasitic loss, but this operator ranks the structures in the same order (Spearman coefficient 0.96).

---

## 3. Methods

### 3.1 Device optics

Dipole emission and stack reflectance are computed with the Chance–Prock–Silbey model [40], extended to birefringent layers [41,42], and with transfer matrices, using measured optical constants and an incoherent substrate [43]. The reference stack is Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / glass ($n$ = 1.51), with a horizontal dipole ratio of 0.865. A second stack replaces Al by Ag, and a low-loss stack follows Ref. [14]. At 550 nm, $A'$ is 0.167, 0.044, and 0.026 for the three stacks (Table S2, Fig. S1). Device optics are computed at each wavelength with dispersive optical constants. The structure matrices are independent of wavelength, because they are obtained by ray optics with a non-dispersive index of 1.51. Spectrally weighted quantities cover 400–700 nm in 10 nm steps and are weighted with the emitter spectrum multiplied by the substrate-coupled fraction $\eta_{\mathrm{sub}}(\lambda)$ (Note S10).

### 3.2 Structures and ray tracing

All structures are index matched to the substrate. The library contains 38 structures, 37 of them non-planar: a planar interface; six non-absorbing Henyey–Greenstein [44] scattering slabs with asymmetry parameter $g$ = 0–0.9 and reduced optical thickness $\tau^*$ = 0.25–4, combined with the exit face by the adding method [45]; nine hexagonal MLAs, namely hemispherical and freeform optima for the total EQE and for four angular bands; 17 random MLAs on a 6 × 6 supercell with independent lenslet profiles (spread $s$), heights, or positions and radii; and five Gaussian-correlated wavy surfaces with rms slopes of 0.2–1.5, which resemble the buckled surfaces used on OLEDs [46]. MLAs and textures are traced with a height-field Monte Carlo ray tracer that reproduces the LightTools return matrix of the hemispherical MLA (escape probability 0.4275 compared with 0.4274). Ray optics is adequate for surface features that are large compared with the wavelength [47]; the lenslets considered here have radii of 10 µm. Numerical parameters and the residual reciprocity asymmetry of the traced matrices change $N_{\mathrm{mix}}$ by at most 0.24, $\xi$ by at most 0.025, and $\eta_{\mathrm{sta}}$ by at most 0.34 pp (Tables S5 and S6, Note S10).

### 3.3 Fabricated devices and finite-geometry model

Devices with two emitters were fabricated on 25.4 mm glass substrates of 0.7 mm thickness. The orange devices are ITO (150 nm) / TAPC (50 nm) / TCTA (10 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / B3PyMPM (45 nm) / LiF (1 nm) / Al (100 nm). The green devices use Ir(ppy)$_2$acac in the same host with nominal transport-layer thicknesses of 105 nm (TAPC) and 42 nm (B3PyMPM). A pixel-defining layer opens a circle of 1 mm radius. For the lens devices, a second 0.7 mm glass carrying a hexagonal array of hemispherical or freeform lenslets (radius 10 µm) on a 3 × 3 mm square centered on the pixel is attached with index-matching oil, and the substrate edges are blackened. The freeform profile differs between the two emitters (Note S13). The EQE and the angular intensity were measured with a photodiode goniometer in 10° steps.

These devices are modeled with the chain extended to lateral position. Between encounters with the front face the light moves laterally by $2t\tan\theta$, where $t$ is the substrate thickness. The back-face reflectance depends on the region of the layout (pixel, Al-covered, or uncovered, each with or without ITO), and light that reaches the blackened edges is lost. The lens kernels are traced with the orientation of the hexagonal lattice resolved, and the goniometer is modeled with its scan plane, rotation axis, and detector size. For an infinite substrate this estimator reproduces $\eta_{\mathrm{sta}}$ to within $2\times10^{-4}$ (Note S13).

The cavity and a single internal efficiency per emitter (0.89 for orange and 0.95 for green) are fixed on the bare devices. For the green stack, reproducing the bare angular profile required transport-layer thicknesses that differ from the nominal ones (65 nm B3PyMPM and 80 nm TAPC; Note S13). The lens devices are then predicted without further parameters. The orientation of the lens lattice relative to the scan plane was not recorded, so predictions are quoted as ranges over lattice orientation, scan plane, rotation axis, and photodiode size.


---

## 4. Results

### 4.1 Where the light is mixed

The effect of repeated passes can be read directly from the angle-resolved series (Al device, 550 nm). On its first arrival at the hemispherical MLA, 28% of the light lies beyond 60° [Fig. 1(e)]. The lenses remove light inside and near the escape cone at every pass and return the rest close to its original angle, so this share grows to 42% after two returns and to 83% after ten, while the maximum of the distribution moves from 63° to 83°. The probability of escape at pass $k$ falls accordingly, from 0.43 at the first encounter to 0.14 at the eleventh, a third of $p$ = 0.42 [Fig. 1(f)]. For the planar interface it drops from 0.35 to 0.05 at the second encounter. Structures that mix within one or two encounters behave differently, but they do not return to $p$ either. The scattering layer settles at 0.316 and the steep wavy surface at 0.366, in both cases within 0.003 of $p_{\mathrm{mix}}$ (0.318 and 0.369) and below $p$ (0.326 and 0.422). A random MLA, with $N_{\mathrm{mix}}$ = 1.9, lies between the two regimes. Light that a structure has mixed in reflection thus escapes with $p_{\mathrm{mix}}$, as Eq. (6) states, and not with the probability $p$ assumed by the fully randomized model.

The size of the penalty varies widely among the structures [Fig. 2(a)]. The ratio $p_{\mathrm{mix}}/p$ is 0.81–0.93 for lens arrays, 0.37–0.90 for wavy surfaces, and 0.14 for the planar interface. Mixing at the device avoids the penalty. The difference in $\eta_{\mathrm{sta}}$ between complete mixing at the device and complete mixing at the escape interface increases with $p - p_{\mathrm{mix}}$, with a correlation coefficient of 0.97 [Fig. 2(b)]. For the planar interface, $\eta_{\mathrm{sta}}$ is 0.384 without mixing, 0.488 with mixing at the interface, 0.770 with mixing at the device, and 0.798 with a lossless transmissive randomizer, the value of the fully randomized model. For lens arrays, whose $\sigma_p^2$ is small, mixing at the device adds 1.5–4.1 pp, compared with at most 1.8 pp at the escape interface. Mixing at the device comes within 2.2 pp of the fully randomized model for all non-planar structures in the Al device and within 1.2 pp in the Ag and low-loss devices, because their reflectance is nearly uncorrelated with $p(\theta)$ ($|p_D - p| < 0.008$, Note S11).

Transmissive mixing does not require an ideal randomizer. A non-absorbing Henyey–Greenstein volume under the planar exit face mixes the light in transmission and also reflects part of it [Fig. 2(d)]. With forward scattering ($g$ = 0.8–0.9) and $\tau^*$ ≈ 0.5, it reaches $\eta_{\mathrm{sta}}$ = 0.764 in the Al device. This is 92% of the way from the planar interface (0.384) to the ideal transmissive value (0.798), and it is far above the 0.488 obtained with complete reflective mixing at the interface. In the Ag and low-loss devices the same volume recovers 97% and 99%. Thicker volumes reflect more of the light and fall back to 0.58 at $\tau^*$ = 4.

The distinction is not specific to OLED substrates. Figure 2(c) shows a slab with an internal isotropic source and an angle-independent back reflectance $R$ = 0.90, for which complete mixing at the back reflector is equivalent to complete transmissive mixing. With a planar exit face and $n$ = 1.5, the extraction efficiency is 0.44 without mixing, 0.60 with complete mixing at the exit face, and 0.87 with complete mixing at the back reflector. At $n$ = 3.5 the three values are 0.08, 0.19, and 0.37. A textured exit face transmits a similar fraction from every internal angle, so its $\sigma_p^2$ is small (Table S28). Hemispherical-MLA and steep wavy faces come within 1.1 pp of the fully randomized result at $n \ge 2.4$, and complete mixing at the exit face changes their efficiency by at most 2.0 pp. A shallow wavy face (rms slope 0.4) stays up to 5.8 pp below the fully randomized result. Mixing at its exit face recovers none of this difference, whereas mixing at the back reflector recovers all of it.

**Fig. 2.** Where the light is mixed. (a) Escape probability of light completely mixed by the escape interface, $p_{\mathrm{mix}}$, against the escape probability $p$ of light with uniform radiance, for all structures. On the dotted line, $p_{\mathrm{mix}} = p$. (b) Difference in $\eta_{\mathrm{sta}}$ between complete mixing at the device and at the escape interface against $p - p_{\mathrm{mix}}$ (Al device, 550 nm). The legend applies to (a), (b), and (e). (c) Extraction efficiency of a slab with an internal isotropic source and an angle-independent back reflectance $R$ against its refractive index, for a planar exit face without mixing, with complete mixing at the exit face, and with complete mixing at the back reflector, which for constant $R$ equals transmissive mixing and the fully randomized model. Symbols show textured exit faces without mixing (Table S28). (d) $\eta_{\mathrm{sta}}$ of non-absorbing Henyey–Greenstein volumes under a planar exit face against the reduced optical thickness $\tau^*$ (Al device, 550 nm), compared with the ideal transmissive randomizer, complete reflective mixing at the planar interface, and the planar interface without mixing. (e) Gain in $\eta_{\mathrm{sta}}$ predicted for complete mixing by the structure at its own $p(\theta)$ (filled symbols) and by the fully randomized model (open symbols) against $N_{\mathrm{mix}}$ (Al device, 550 nm). The planar interface (10.4 and 41.4 pp) is off scale. *(`fig2_complete_mixing.png`)*

### 4.2 Gain from complete mixing by the structure

How much would each structure gain if it mixed the light completely while keeping its own $p(\theta)$? At 550 nm, Eq. (7) changes $\eta_{\mathrm{sta}}$ of the 37 non-planar structures by −0.4 to +1.8 pp in the Al device. The largest increase is 2.0 pp in the Ag device and 0.4 pp in the low-loss device [Fig. 2(e)]. Weighted with the emission spectrum (peak at 572 nm), the largest increases are 1.1, 1.4, and 2.2 pp for the Al, Ag, and low-loss devices. The fully randomized model predicts more, by 0.6–2.9 pp for lens arrays, by 17 pp for the shallowest wavy surface, and by 31 pp for the planar interface. Mixing by the device, in comparison, adds 1.1–3.3 pp for lens arrays in the Al device under the same spectral weighting (Table S29).

The gain is larger when the cavity emits into angles that the structure transmits poorly. At 500 nm the first-pass escape probability of the hemispherical MLA is 0.36 in the Al device and 0.33 in the Ag device, well below $p$ = 0.42. Complete mixing would then raise $\eta_{\mathrm{sta}}$ by 5.7–6.4 pp (Al) and 8.4–9.4 pp (Ag), and random lenses and wavy surfaces realize most of this increase (Table S10). Reflective mixing is useful, in other words, when the first-pass angular distribution does not match $p(\theta)$.

Two wavy surfaces (rms slopes 0.4 and 0.7) exceed $\eta_{\mathrm{mix}}$, by 0.42 and 0.01 pp. The second value is below the numerical uncertainty of 0.34 pp. These are the only structures in which the escape probabilities at successive encounters are anticorrelated (Note S10).

### 4.3 Accuracy and limits of the reduced model

We compared four models with the full series of Eq. (2) for 342 cases, comprising 38 structures, three devices, and three wavelengths. The fully randomized model overestimates $\eta_{\mathrm{sta}}$ of the planar interface by up to 60 pp (Ag device, 500 nm) and that of the hemispherical MLA by 2–3 pp [Fig. S2(a)]. The opposite assumption, that every ray keeps its angle, underestimates the scattering layers [Fig. S2(b)]. A reduced model with $m = 1-\lambda_2$ and isotropic redistribution fits the lens arrays but overestimates shallow wavy surfaces by up to 10.5 pp [Fig. S2(c)]. These surfaces have a small $\lambda_2$ but return trapped light mostly to other trapped angles, which reciprocal redistribution captures without a further parameter. With it, the mean absolute error falls to 0.30 pp (mean signed error −0.24 pp) and the largest error to 3.0 pp [Fig. 3(a)]. The error is not concentrated in any device or wavelength [Fig. 3(b), Table S4].

Direct ray tracing in $(\theta, \phi)$ agrees with the full series to within 0.6 pp for eight structures (Note S9, Table S19). The remaining differences vanish when the azimuth is randomized at each pass, which identifies them with the azimuthal correlation that the chain neglects. For nine test structures that were not used in developing the model (spherical caps, paraboloids, Gaussian bumps, a cone–hemisphere mixture, wavy surfaces with a ring spectrum, and two cone arrays), the mean error is 0.48 pp and the largest 3.3 pp for seven of them. For the two cone arrays the error reaches 20 pp.

The two failures point to conditions under which a single rate $1-\lambda_2$ cannot describe the circulating light [Fig. 3(c)]. Both can be checked from $M$ and the first-pass distribution without evaluating the series. The chain must be irreducible. In the 45° cone array, light within 0–4° (0.5% of the étendue) is retroreflected into the same range, a closed set that gives $\lambda_2 = 1$ although all other angles mix. In addition, the slowest mode must carry a significant share of the returned light, which we quantify by the relaxation rate weighted with the returned light, $\lambda_{\mathrm{eff}}$ (Note S9). The 60° cones have $\lambda_2$ = 0.985 but $\lambda_{\mathrm{eff}}$ = 0.43. All 45 other structures have $\lambda_2 - \lambda_{\mathrm{eff}} \le 0.35$ and a largest error of 3.3 pp.

Because these conditions were formulated after the out-of-sample test, we tested them prospectively on twelve further structures, namely cones with side angles of 30–70°, square and inverted pyramids, and one-dimensional prism gratings, and recorded the predictions before evaluating the series (Note S9). The irreducibility condition held. The 45° pyramids and the 45° prisms form closed sets, were predicted to fail, and have errors of 21 and 23 pp. The gap condition did not discriminate. It flagged two structures that the model describes to within 2.5 pp, and the eight structures that passed both conditions have a mean error of 1.4 pp but a largest error of 7.5 pp, for 40° cones. The reduced model is therefore less accurate for faceted structures than for the curved lenses, scatterers, and wavy surfaces of the library, although its mean error stays below that of the fully randomized model (2.8 pp). Direct traces confirm the full series for the 60° cones to within 0.6 pp, so these failures belong to the one-rate reduction and not to the angle-resolved chain.

**Fig. 3.** Accuracy of the reciprocal reduced model. (a) $\eta_{\mathrm{sta}}$ from the reduced model, with $m = 1-\lambda_2$ and $q \propto w(1-p)$, against the full series for 342 cases. (b) Error in (a) for each device (LL, low loss) and wavelength. Points are offset horizontally for visibility, and bars give the mean absolute error. (c) Largest error over the nine device–wavelength cases against $\lambda_2 - \lambda_{\mathrm{eff}}$ (Al device, 550 nm) for the 38 library structures, nine out-of-sample structures, and twelve structures tested prospectively. Crosses mark reducible chains (45° cones, pyramids, and prisms), and the shaded region marks the gap threshold identified after the out-of-sample test. *(`fig3_reduced_model.png`)*

### 4.4 Escape probability and mixing number of the structures

The structure classes occupy separate regions of the $(p, N_{\mathrm{mix}})$ plane [Fig. 4(a)]. The planar interface has $p$ = 0.40 and $N_{\mathrm{mix}} = \infty$. Scattering layers relax within 1.1–2.2 encounters but have lower escape probabilities, $p$ = 0.19–0.39. Lens arrays have the highest escape probabilities, $p$ = 0.42–0.43, with $N_{\mathrm{mix}}$ = 6–8 for the hemispherical and 3–7.5 for the freeform arrays.

Randomizing a lens array lowers $N_{\mathrm{mix}}$ without lowering $p$ [Fig. 4(b), Table S3]. Independent lenslet profiles reduce $N_{\mathrm{mix}}$ from 7.0 to 1.9, a position jitter of 30% of the pitch reduces it to 3.5, and a height spread of ±60% only to 4.1 (Table S5). Wavy surfaces with rms slopes of 0.4–1.5 reach $N_{\mathrm{mix}}$ = 1.1–1.2, as low as a scattering layer.

The late-pass escape ratio carries information that $N_{\mathrm{mix}}$ does not [Fig. 4(c)]. For lens arrays, $\xi$ decreases from 0.80 to 0.39 as $N_{\mathrm{mix}}$ increases from 2 to 8, because light that keeps its angle accumulates where the lenses transmit least. The shallowest wavy surface ($\xi$ = 0.35) and a thin isotropic scattering layer ($\xi$ = 0.54) lie below this trend at $N_{\mathrm{mix}} \approx 2$. They relax quickly but return trapped light mainly to other trapped angles, and they show the largest residual errors in Fig. 3(a).

The contours of $\eta_{\mathrm{sta}}$ in Fig. 4(a) are nearly vertical, so the efficiency is set mainly by $p$. For the $p(\theta)$ shapes of lenses, scattering layers, and wavy surfaces, the full range of $N_{\mathrm{mix}}$ changes $\eta_{\mathrm{sta}}$ only as much as a change of $p$ by 0.02–0.05 (Table S8). Randomizing the lenslet profiles, for example, reduces $N_{\mathrm{mix}}$ from 7.0 to 1.9 but changes $\eta_{\mathrm{sta}}$ by less than 1 pp in the Al device at 550 nm, because the gain from faster mixing (+0.9 pp) is offset by the accompanying change of $p(\theta)$ (−2.2 pp, Table S9).

**Fig. 4.** Escape probability and mixing number of the structures. (a) Escape probability and angular mixing number of all non-planar structures, with contours of $\eta_{\mathrm{sta}}$ (reduced model, Al device, 550 nm, hemispherical $p(\theta)$ scaled to each $p$). The hemispherical and freeform MLAs are the optima of the benchmark (Table S15). (b) $N_{\mathrm{mix}}$ against each randomization parameter, normalized to its maximum. Dotted lines mark the hemispherical MLA and a scattering layer. (c) $\xi$ against $N_{\mathrm{mix}}$. The planar interface ($N_{\mathrm{mix}} = \infty$, $\xi$ = 0) is drawn at the right edge. *(`fig4_structure_map.png`)*

### 4.5 Fabricated devices

The lens patch of the fabricated devices is only 3 mm wide [Fig. 5(a)], on a substrate that is 1.4 mm thick once the lens glass is attached. Most light returned by the lenses lands outside the patch, where the planar face traps it until it reaches the blackened edges. The lenses therefore act mainly on first-pass light, and the experiment tests how the model redistributes that light in angle.

The model predicts the angular signature of both lens types for both emitters without adjustment [Figs. 5(b)–5(d), Table 2]. The hemispherical array concentrates the emission toward the normal. Its angular factor, defined as the EQE divided by the EQE of a Lambertian emitter with the same on-axis intensity, is predicted as 0.62–0.68 for the orange and 0.65–0.71 for the green devices and measured as 0.66 ± 0.04 and 0.71 ± 0.02. The freeform arrays spread the first-pass light, and the bare devices are close to Lambertian. All device means lie inside the predicted ranges except that of the orange freeform devices, which is 0.03 above its range. The mean angular profiles deviate from the center of the predicted bands by at most 0.04 rms between 10° and 80° [Figs. 5(b) and 5(c)]. Values for the individual devices are listed in Tables S32 and S33.

**Table 2.** Angular factor of the fabricated devices. The model range covers lattice orientation, scan plane, rotation axis, and photodiode size for the 3 × 3 mm lens patch. Measured values are the mean ± standard deviation over devices, or over scans where only one device was available.

| Emitter | Device | Model | Measured | Sample |
|---|---|---|---|---|
| Orange | Bare | 0.99–1.02 | 1.01 ± 0.05 | 5 devices |
| Orange | Hemispherical MLA | 0.62–0.68 | 0.66 ± 0.04 | 5 devices |
| Orange | Freeform MLA | 0.92–1.04 | 1.07 ± 0.07 | 4 devices |
| Green | Bare | 1.00–1.02 | 1.02 | 2 devices |
| Green | Hemispherical MLA | 0.65–0.71 | 0.71 ± 0.02 | 1 device, 6 scans |
| Green | Freeform MLA | 0.83–0.91 | 0.87 ± 0.04 | 3 devices |

The measurements also show that on-axis metrics overstate what a small lens patch achieves [Fig. 5(e)]. The hemispherical devices have 1.22 (orange) and 1.21 (green) times the on-axis EQE of the bare devices but only 0.81 and 0.85 times their total EQE. The lossless model predicts the same contrast at a higher level, 1.39–1.55 on axis against 0.94–1.01 in total, because the absolute EQE of all lens devices is 84–91% of the lossless prediction. The role of the lateral extent of an outcoupling structure, and of the area over which light is collected, has been pointed out for laboratory-size samples [22] and for pixelated displays [11].

A single angle-independent loss in the attached lens glass accounts for this offset. A transmission of 0.925 per crossing of the oil-coupled interface, fitted to the total EQE of the orange hemispherical devices and the only parameter fitted to lens devices, predicts 0.89, 0.87, and 0.90 of the lossless EQE for the orange freeform, green hemispherical, and green freeform devices, compared with 0.90, 0.84, and 0.91 measured, and changes the angular factors by at most 0.02 (Note S13). With this loss, the predicted on-axis and total ratios are 1.19–1.33 and 0.81–0.87 for the hemispherical devices and 0.83–1.03 and 0.84 for the freeform devices, against measured values of 1.21–1.22 and 0.81–0.85, and 0.80–1.00 and 0.85 [filled boxes in Fig. 5(e)]. Variations of the lens shape cannot explain the offset, because those that lower the predicted on-axis gain to the measured value raise the angular factor above its measured range (Table S35). For a 10 × 10 mm patch the model predicts 1.28–1.77 times the bare EQE.

**Fig. 5.** Fabricated devices and finite-geometry model. (a) Substrate layout with the pixel, the 3 × 3 mm lens patch, and the 10 × 10 mm patch that is used only in the model. (b), (c) Normalized angular intensity of the orange and green devices. Bands show the model range over the measurement geometries for the 3 × 3 mm patch, with the band center as a line. Symbols are measurements: in (b), one bare device, the mean and range of two freeform devices, and the mean ± standard deviation of five hemispherical devices; in (c), the mean ± standard deviation of the tilt-corrected scans of the bare (seven scans, two devices), hemispherical (six scans, one device), and freeform (twelve scans, three devices) devices. (d) Angular factor: measured mean ± standard deviation (Table 2), model range for the 3 × 3 mm patch, and model value for the 10 × 10 mm patch. (e) Ratio of lens to bare devices for the Lambertian-equivalent on-axis EQE (axial) and for the total EQE. Dashed outlines show the lossless model, and filled boxes the model with a transmission of 0.925 per crossing of the lens-glass interface, fitted to the total EQE of the orange hemispherical devices. Measured values are the mean ± standard deviation over devices. *(`fig5_experiment.png`)*

### 4.6 Robustness

Table 3 collects the conditions under which the central result was tested. In every case, reflective mixing by the escape interface adds little, and mixing in transmission or at an angle-independent reflector recovers the fully randomized model. The exception already noted is a first-pass distribution that does not match $p(\theta)$, as at 500 nm.

**Table 3.** Robustness of the central result. Gains are changes of $\eta_{\mathrm{sta}}$ in pp for complete mixing, for the 37 non-planar structures unless stated otherwise.

| Condition | Result | Source |
|---|---|---|
| Al, Ag, and low-loss devices, 550 nm | Reflective mixing by the structure: −0.4 to +1.8 (Al), at most 2.0 (Ag) and 0.4 (low loss). Mixing at the device: +1.5 to +4.1 for lens arrays (Al) | Secs. 4.1 and 4.2 |
| Spectral weighting, 400–700 nm | Reflective mixing: at most 1.1, 1.4, and 2.2 (Al, Ag, low loss). Mixing at the device: 1.1–3.3 for lens arrays (Al) | Table S29 |
| First pass not matched to $p(\theta)$, 500 nm | Reflective mixing: 5.7–6.4 (Al) and 8.4–9.4 (Ag) for the hemispherical MLA | Table S10 |
| Exit medium with $n$ = 1.5–3.5 (slab, $R$ = 0.90) | Textured faces within 1.1 of the fully randomized model at $n \ge 2.4$; mixing at the exit face at most 2.0 | Sec. 4.1, Table S28 |
| Finite substrate (25 mm, absorbing edges) | Reflective mixing: −5.2 to +0.8. Transmissive mixing: −1.0 to +22.7 (six structures, three devices) | Fig. S10 |
| Partly diffuse electrode (diffuse fraction 0–1) | Gain nearly linear in the diffuse fraction, at most 2.3 | Table S30 |
| Numerical parameters and reciprocity residual | $\eta_{\mathrm{sta}}$ changes by at most 0.34, $N_{\mathrm{mix}}$ by at most 0.24 | Tables S5 and S6 |
| Separate s and p channels | At most 0.4 | Table S7 |
| Direct $(\theta,\phi)$ ray tracing | Within 0.6 of the full series | Note S9 |
| Angular cutoff of the spectral quantities, 80–90° | Mean error of the reduced model changes by less than 0.3 | Note S1 |

---

## 5. Discussion

### 5.1 Reflective and transmissive mixing

Reciprocity couples angular redistribution to escape. A surface that mixes light in reflection sends it back in proportion to what each angle reflects, so the mixed light is concentrated where escape is least likely, and its escape probability falls below $p$ by the variance term of Eq. (6). The fully randomized model assumes the opposite weighting, which only a transmissive mixer provides. This explains why faster mixing by a refractive structure does not by itself raise the extraction. In the structures examined here it usually comes with a reshaped or lower $p(\theta)$, as for the randomized lens arrays (Table S9) and for scattering layers, whose efficiency peaks at an intermediate scattering strength (Fig. 2(d); see also Ref. [18]). A volume with little backscattering can approach the fully randomized model. A refractive structure that mixes through its own escape interface cannot.

### 5.2 Limits of lens-shape optimization

A hemispherical array already reaches an escape probability close to the étendue limit for uniform radiance, $1/n^2$. Changing the lens profile then mainly moves $N_{\mathrm{mix}}$ and the shape of $p(\theta)$, directions along which the efficiency contours of Fig. 4(a) are flat. The anisotropic emission of the cavity would in principle allow more. For the Al stack at 550 nm, matching the étendue of the escape cone to the brightest internal angles gives a first-pass bound of 0.540, compared with $1/n^2$ = 0.439 (Fig. S4, Note S4); bounds of this kind follow from étendue conservation [48,49]. The lenses nevertheless transmit almost the same power as for uniform radiance (0.43), because they fill each outgoing direction from a broad range of internal angles. In a controlled LightTools benchmark on the device platform of Refs. [11,12], band-optimized freeform lenses exceed an equally optimized hemisphere by at most 7% in one angular band and by 1.3% in total EQE (Note S12). Approaching the bound appears to require a surface whose transmission is angle selective. The angular composition of the output is also nearly independent of the structure. The fraction emitted into 40–60° is 0.32–0.36, close to the Lambertian value of 0.337 (Fig. S3, Table S11).

### 5.3 Implications for design

Before a structure is optimized for mixing, Eq. (6) gives an estimate of what mixing can add. For the lens arrays studied, $\sigma_p^2/(1-p)$ is 0.03–0.08 and the gain is at most about 2 pp. Larger gains require that the mixing be placed where reciprocity does not penalize it. One option is a weakly absorbing, forward-scattering volume. With $g \ge 0.8$ and $\tau^* \approx 0.5$ it recovers 92–99% of the transmissive gain over a planar interface. The other is the device itself, provided that its reflectance is uncorrelated with $p(\theta)$, as it is for the three stacks considered here; for lens arrays in the Al device this adds 1.5–4.1 pp.

When structures are compared, $p$ should be ranked first. $N_{\mathrm{mix}}$ and $\xi$ matter mainly when the first-pass distribution of the cavity does not match $p(\theta)$, where mixing can add 5–9 pp. Random MLAs and wavy surfaces with an rms slope of one or more combine a high $p$ with a low $N_{\mathrm{mix}}$ and a high $\xi$. Their spectrally weighted advantage over periodic MLAs is nevertheless only 0.2–0.5 pp in the Al and Ag devices, comparable to the numerical uncertainty, and 2 pp in the low-loss device. Irregular microlens arrays remain attractive for another reason, since they suppress the angular color shift of microcavity devices [50]. The $(p, N_{\mathrm{mix}})$ plane can also serve as a diagnostic for numerical optimizers, including the learning-based ones [9,10], because it shows whether a search is moving along or across the efficiency contours.

The one-rate reduced model should be used only for irreducible chains, and errors of several percentage points should be expected for faceted structures (Sec. 4.3). The full matrix costs little more and has neither restriction.

The fabricated devices add two practical points. A lens patch must extend several substrate thicknesses beyond the pixel to extract recycled light. For the 1.4 mm substrate of Sec. 4.5 the model predicts that a hemispherical patch of 3 mm lowers the EQE, that about 6 mm gives half the gain of a 20 mm patch, and that 10 mm gives 80% of it (Table S36). The angle-integrated EQE should also be reported, because in the first-pass regime the on-axis EQE overstates the effect of a lens patch by a factor of about 1.5.

### 5.4 Scope and limitations

The analysis applies to passive external structures that are index matched to the substrate and large enough for ray optics. It does not cover the extraction of waveguided modes by internal structures or the change of the cavity that such structures cause, nor angle-selective filters, absorbing scatterers, or wavelength converters. Except in Secs. 2.5 and 4.1, the device is treated as an angle-preserving reflector, as in earlier substrate models [2,18]. For a partly diffuse electrode the gain grows nearly linearly with the diffuse fraction (Table S30), but the product of the device kernel with $B_R$ need not have a real spectrum. The experiment probes the first-pass regime of small lens patches (Sec. 4.5); the multi-pass gains of Secs. 4.1–4.4 rest on the calculations. Light trapping in solar cells is the reciprocal problem, in which a low escape probability is the aim, and is not analyzed here.

---

## 6. Conclusion

Optical reciprocity ties the value of angular mixing to its placement. Light mixed in reflection by the escape interface escapes with probability $p - \sigma_p^2/(1-p)$, light mixed in transmission through a lossless volume with $p$, and light mixed at the device with $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$. The fully randomized model of statistical ray optics corresponds to the second case. For 37 scattering, microlens, random-lens, and wavy structures, complete mixing by the structure raises the spectrally weighted extraction by at most about 2 pp, whereas the fully randomized model predicts up to 19 pp at 550 nm, and a forward-scattering volume recovers 92–99% of the transmissive gain over a planar interface. These statements hold for three device stacks, for exit media with refractive indices up to 3.5, and for finite substrates. The reciprocal reduced model predicts multi-pass extraction with a mean absolute error of 0.30 pp for curved and random structures and fails, as predicted, for reducible chains. On fabricated devices with two emitters, the finite-geometry model reproduces the angular factors to within 0.03 of the device means without adjustment and shows that small lens patches raise the on-axis EQE while lowering the total extraction. Substantially higher extraction will require mixing in transmission or at the device, lens patches much wider than the substrate thickness, or angle-selective structures.

---

## Funding

[TBD]

## Acknowledgment

[TBD]

## Disclosures

The authors declare no conflicts of interest.

## Data availability

Data underlying the results presented in this paper are not publicly available at this time but may be obtained from the authors upon reasonable request.

## Supplemental document

See Supplement 1 for supporting content.

---

## References

1. W. Brütting, J. Frischeisen, T. D. Schmidt, et al., "Device efficiency of organic light-emitting diodes: progress by improved light outcoupling," Phys. Status Solidi A **210**, 44–65 (2013).

2. H. Greiner, "Light extraction from organic light emitting diode substrates: simulation and experiment," Jpn. J. Appl. Phys. **46**, 4125–4137 (2007).

3. S. Möller and S. R. Forrest, "Improved light out-coupling in organic light emitting diodes employing ordered microlens arrays," J. Appl. Phys. **91**, 3324–3327 (2002).

4. Y. Qu, J. Kim, C. Coburn, et al., "Efficient, nonintrusive outcoupling in organic light emitting devices using embedded microlens arrays," ACS Photonics **5**, 2453–2458 (2018).

5. J. Song, K.-H. Kim, E. Kim, et al., "Lensfree OLEDs with over 50% external quantum efficiency via external scattering and horizontally oriented emitters," Nat. Commun. **9**, 3207 (2018).

6. J. W. Levell, S. Harkema, R. K. Pendyala, et al., "Device reflectivity as a simple rule for predicting the suitability of scattering foils for improved OLED light extraction," Proc. SPIE **8829**, 88291L (2013).

7. M. Kovačič, "Design of external microtextures for efficient light outcoupling in OLEDs with different preferential orientation of emission dipoles," Opto-Electron. Rev. **30**, e141542 (2022).

8. Y. Sun and S. R. Forrest, "Organic light emitting devices with enhanced outcoupling via microlenses fabricated by imprint lithography," J. Appl. Phys. **100**, 073106 (2006).

9. S. Kim, J. M. Shin, J. Lee, et al., "Inverse design of organic light-emitting diode structure based on deep neural networks," Nanophotonics **10**, 4533–4541 (2021).

10. Z.-R. Zhao, K.-P. Song, H. Liu, et al., "Deep Q-network-guided TiO$_2$/NOA61 microstructures for light outcoupling and angular emission redistribution in blue top-emitting OLEDs," Adv. Funct. Mater. **36**, e78939 (2026).

11. M. Kim, J. Kim, and S. Yoo, "Near-planar light outcoupling structures with finite lateral dimensions for ultra-efficient and optical crosstalk-free OLED displays," Nat. Commun. **16**, 11606 (2025).

12. J. Kim, E. Kim, J. Park, et al., "Toward near-foldable surface light sources with ultimate efficiency: ultrathin substrates embedded with micron-scale inverted lens arrays," ACS Photonics **10**, 1775–1782 (2023).

13. W. Li, Y.-Q. Li, Y. Shen, et al., "Releasing the trapped light for efficient silver nanowires-based white flexible organic light-emitting diodes," Adv. Opt. Mater. **7**, 1900985 (2019).

14. J. Kim, S. Lee, S. Yoo, et al., "Minimizing round-trip loss for near-unity quantum efficiency in organic light-emitting diodes," arXiv, arXiv:2610.12065 (2026).

15. E. Yablonovitch, "Statistical ray optics," J. Opt. Soc. Am. **72**, 899–907 (1982).

16. I. Schnitzer, E. Yablonovitch, C. Caneau, et al., "30% external quantum efficiency from surface textured, thin-film light-emitting diodes," Appl. Phys. Lett. **63**, 2174–2176 (1993).

17. M. A. Green, "Lambertian light trapping in textured solar cells and light-emitting diodes: analytical solutions," Prog. Photovoltaics **10**, 235–241 (2002).

18. J. J. Shiang and A. R. Duggal, "Application of radiative transport theory to light extraction from organic light emitting diodes," J. Appl. Phys. **95**, 2880–2888 (2004).

19. J. Eisenlohr, N. Tucher, O. Höhn, et al., "Matrix formalism for light propagation and absorption in thick textured optical sheets," Opt. Express **23**, A502–A518 (2015).

20. A. David, "Surface-roughened light-emitting diodes: an accurate model," J. Disp. Technol. **9**, 301–316 (2013).

21. S. Altazin, C. Reynaud, U. M. Mayer, et al., "Simulations, measurements, and optimization of OLEDs with scattering layer," SID Symp. Dig. Tech. Pap. **46**, 564–567 (2015).

22. M. Kovačič, D. Samigullina, F. Bouchard, et al., "Analysis and optimization of light outcoupling in OLEDs with external hierarchical textures," Opt. Express **29**, 23701 (2021).

23. H. W. Deckman, C. B. Roxlo, and E. Yablonovitch, "Maximum statistical increase of optical absorption in textured semiconductor films," Opt. Lett. **8**, 491–493 (1983).

24. S. Manzoor, M. Filipič, A. Onno, et al., "Visualizing light trapping within textured silicon solar cells," J. Appl. Phys. **127**, 063104 (2020).

25. C. Battaglia, M. Boccard, F.-J. Haug, et al., "Light trapping in solar cells: when does a Lambertian scatterer scatter Lambertianly?" J. Appl. Phys. **112**, 094504 (2012).

26. M. C. W. van Rossum and T. M. Nieuwenhuizen, "Multiple scattering of classical waves: microscopy, mesoscopy, and diffusion," Rev. Mod. Phys. **71**, 313–371 (1999).

27. D. A. Levin and Y. Peres, *Markov Chains and Mixing Times*, 2nd ed. (American Mathematical Society, 2017).

28. J. N. Darroch and E. Seneta, "On quasi-stationary distributions in absorbing discrete-time finite Markov chains," J. Appl. Probab. **2**, 88–100 (1965).

29. E. G. Altmann, J. S. E. Portela, and T. Tél, "Leaking chaotic systems," Rev. Mod. Phys. **85**, 869–918 (2013).

30. J.-W. Ryu, S.-Y. Lee, C.-M. Kim, et al., "Survival probability time distribution in dielectric cavities," Phys. Rev. E **73**, 036207 (2006).

31. K. Clauß, E. G. Altmann, A. Bäcker, et al., "Structure of resonance eigenfunctions for chaotic systems with partial escape," Phys. Rev. E **100**, 052205 (2019).

32. E. Veach, "Robust Monte Carlo methods for light transport simulation," Ph.D. thesis, Stanford University (1997).

33. C. Soler and K. Subr, "Spectral theory of light transport operators," ACM Trans. Graph. **45**, 1–18 (2026).

34. J. J. Koenderink and A. J. van Doorn, "Geometrical modes as a general method to treat diffuse interreflections in radiometry," J. Opt. Soc. Am. **73**, 843 (1983).

35. F. J. J. Clarke and D. J. Parry, "Helmholtz reciprocity: its validity and application to reflectometry," Light. Res. Technol. **17**, 1–11 (1985).

36. M. Ashikhmin and P. Shirley, "An anisotropic Phong BRDF model," J. Graph. Tools **5**, 25–32 (2000).

37. C. Kelemen and L. Szirmay-Kalos, "A microfacet based coupled specular-matte BRDF model with importance sampling," in *Eurographics 2001 – Short Presentations* (Eurographics Association, 2001).

38. K. Neyts and A. Ullan Nieto, "Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices," J. Opt. Soc. Am. A **23**, 1201 (2006).

39. R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed. (Cambridge University Press, 2012).

40. R. R. Chance, A. Prock, and R. Silbey, "Molecular fluorescence and energy transfer near interfaces," Adv. Chem. Phys. **37**, 1–65 (1978).

41. J. A. E. Wasey, A. Safonov, I. D. W. Samuel, et al., "Effects of dipole orientation and birefringence on the optical emission from thin films," Opt. Commun. **183**, 109–121 (2000).

42. C.-K. Moon, S.-Y. Kim, J.-H. Lee, et al., "Luminescence from oriented emitting dipoles in a birefringent medium," Opt. Express **23**, A279 (2015).

43. S. Mladenovski, S. Hofmann, S. Reineke, et al., "Integrated optical model for organic light-emitting devices," J. Appl. Phys. **109**, 083114 (2011).

44. L. G. Henyey and J. L. Greenstein, "Diffuse radiation in the galaxy," Astrophys. J. **93**, 70–83 (1941).

45. H. C. van de Hulst, *Multiple Light Scattering: Tables, Formulas, and Applications* (Academic Press, 1980).

46. W. H. Koo, S. M. Jeong, F. Araoka, et al., "Light extraction from organic light-emitting diodes enhanced by spontaneously formed buckles," Nat. Photonics **4**, 222–226 (2010).

47. K. Tang, R. A. Dimenna, and R. O. Buckius, "Regions of validity of the geometric optics approximation for angular scattering from very rough surfaces," Int. J. Heat Mass Transfer **40**, 49–59 (1997).

48. J. Muschaweck and C. Wiesmann, "Extracting light out of LEDs," Adv. Opt. Technol. **2**, 291–300 (2013).

49. A. Lenef, A. Piquette, and J. Kelso, "Thermodynamics of light extraction from luminescent materials," ECS J. Solid State Sci. Technol. **7**, R3211–R3226 (2018).

50. J. Lim, S. S. Oh, D. Y. Kim, et al., "Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens array," Opt. Express **14**, 6564 (2006).
