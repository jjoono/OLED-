# Reflective and transmissive angular mixing affect substrate light extraction differently in organic light-emitting diodes

**Authors** — [TBD]

---

## Abstract

Angular mixing is often assumed to help extract the light trapped in the substrate of organic light-emitting diodes, as in the fully randomized model of statistical ray optics. Using optical reciprocity and angle-resolved transport calculations, we show that the benefit depends on where the light is mixed. An escape interface that mixes light in reflection sends it back mostly to angles from which it cannot escape, so the mixed light escapes less often, by an amount set by the angular variance of the escape probability. Mixing in transmission through a lossless volume has no such penalty and reproduces the fully randomized model. For 37 microlens, scattering, and textured structures, complete mixing by the structure raises the spectrally weighted extraction by at most 2.2 percentage points, unless the cavity emits mainly into angles that the structure transmits poorly. A forward-scattering volume recovers 92–99% of the fully randomized gain over a planar interface. A one-parameter model predicts 342 full calculations with a mean error of 0.30 percentage points but is less accurate for faceted structures. Microlens devices with two emitters confirm the predicted angular redistribution of first-pass light.

---

## 1. Introduction

A planar organic light-emitting diode (OLED) emits only part of the light it generates. The rest is lost to waveguided and surface-plasmon modes, to absorption in the electrodes, and to total internal reflection at the substrate–air interface [@brutting2013;@greiner2007]. The last loss can be reduced without changing the electrical stack, by adding microlens arrays (MLAs) [@moller2002;@qu2018], scattering layers [@song2018;@levell2013], or textured films [@kovacic2022] to the substrate. How extraction depends on lens profile, aspect ratio, and fill factor is well known [@moller2002;@sun2006], and deep learning is now used to design both the thin-film stack [@kim2021] and the outcoupling structure [@zhao2026]. Even so, the external quantum efficiencies (EQEs) reported with such structures mostly lie between about 50% and 60% [@song2018;@kim2025;@kim2023;@li2019]. A recent analysis traces this limit to the loss that occurs each time unextracted light returns to the device [@kim2026].

Two kinds of model describe this multi-pass process. Statistical ray optics assumes that the light in the substrate is fully randomized at every pass [@yablonovitch1982;@schnitzer1993;@green2002]. The substrate-to-air efficiency is then $\eta = p/[p+(1-p)A']$, where $p$ is the escape probability of light with uniform radiance and $A'$ is the round-trip loss [@kim2026]. Angle-resolved models instead follow the angular power distribution from pass to pass with redistribution matrices from ray tracing or radiative transfer [@greiner2007;@kim2026;@shiang2004;@eisenlohr2015;@david2013;@altazin2015;@kovacic2021]. They are accurate, but they give one number per structure and no rule for comparing structures. How far a real texture is from the randomized limit has mostly been studied case by case: through the absorption enhancement of textured films [@deckman1983], pass-resolved ray tracing of textured silicon [@manzoor2020], and the conditions under which a rough interface scatters like a Lambertian surface [@battaglia2012].

The two descriptions differ most where it matters for design. A hemispherical MLA returns most unextracted light close to the angle at which it arrived [@kim2026]. Light that survives several passes therefore collects at angles that the lenses transmit poorly, and the fully randomized model overestimates extraction by an amount that varies from structure to structure. This suggests designing structures that mix angles faster. Whether faster mixing actually helps, and by how much, has not been shown.

Here we show that the answer depends on where the light is mixed, and that optical reciprocity sets this dependence. For a passive structure, reciprocity makes the return matrix symmetric once each angle is weighted by its étendue. An escape interface that mixes light in reflection must then send light back to each angle in proportion to how much that angle reflects. Completely mixed light therefore escapes with probability $p - \sigma_p^2/(1-p)$, where $\sigma_p^2$ is the angular variance of the escape probability. Mixing in transmission through a lossless volume has no such penalty. Mixing at the device is penalized only if the device reflectance is correlated with the escape probability. The fully randomized model is thus a model of transmissive mixing. It does not describe the reflective mixing that refractive outcoupling structures provide.

Section 2 derives these results and a reduced model with one mixing parameter. This parameter is set by the second eigenvalue of the return matrix and defines an angular mixing number. The mixing number plays the role, per pass, of the transport factor $1/(1-g)$ in radiative transfer [@vanrossum1999] and of the relaxation time of a reversible Markov chain [@levin2017]. The largest eigenvalue gives the escape probability of light that has circulated many times. This is the quasi-stationary distribution of an absorbing chain [@darroch1965] and corresponds to the escape rate of an open dynamical system [@altmann2013;@ryu2006;@clauss2019]. Such spectral descriptions are established for light-transport operators [@veach1997;@soler2026;@koenderink1983], but their consequences for light extraction have not been worked out. Section 4 evaluates 38 structures in three device stacks. It tests the reduced model against full transport, direct ray tracing, and structures chosen after the model was fixed, and it compares a finite-geometry version of the model with fabricated MLA devices of two emitters. Section 5 discusses the étendue limit and what the results mean for design.

---

## 2. Theory

### 2.1 Multi-pass extraction as a Markov chain on angle

We divide the polar angle inside the substrate into 1° bins $k$, from 0° to 90°. The étendue of bin $k$ is $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$, so uniform radiance means power proportional to $w_k$. Each time light meets the structure, a fraction $p_k$ of the light in bin $k$ leaves into air, and a fraction $B_R(i\to r)$ of the light in bin $i$ returns into bin $r$ [Fig. 1(a)]. The structures are non-absorbing, so on the full 0–90° grid every incident bin satisfies $p_i + \sum_r B_R(i\to r) = 1$. We write the return matrix as

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}), \tag{1}
$$

which separates how much light returns from where it goes. Each column of $M$ sums to one, so $M$ is the transition matrix of a Markov chain on angle. The étendue-weighted mean $p = \sum_k w_k p_k/\sum_k w_k$ is the escape probability of light with uniform radiance.

The light from the cavity enters through the launch vector $\mathbf{v}^{(0)}$. Its entry $v^{(0)}_k$ is the power that the emitter couples into the substrate in bin $k$, summed over azimuth and over s and p polarization. It is computed on a fine angular grid and summed into the 1° bins, and it is normalized so that $\sum_k v^{(0)}_k = 1$. All efficiencies below are therefore fractions of the power coupled into the substrate. The device reflects returned light with its angle-resolved reflectance $R(\theta)$ and, unless stated otherwise, keeps its angle. If $\mathbf{v}^{(j)}$ is the light that meets the structure for the $(j+1)$th time, then

$$
\mathbf{v}^{(j+1)} = \mathrm{diag}\,\mathbf{R}\; B_R\,\mathbf{v}^{(j)} , \qquad \eta_{\mathrm{sta}} = \sum_{j\ge 0}\mathbf{p}^{\mathsf T}\mathbf{v}^{(j)} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} . \tag{2}
$$

Here $j = 0$ is the first encounter, and $\eta_{\mathrm{sta}}$ is the substrate-to-air extraction efficiency [@kim2026]. The escape probability at encounter $j$ is $p_j = \mathbf{p}^{\mathsf T}\mathbf{v}^{(j)}/\sum_k v^{(j)}_k$. All losses enter through $R(\theta)$. The scalar $A' = 1-\sum_k w_k R_k/\sum_k w_k$ is used only to compare with the fully randomized model.

The model is scalar and averages over azimuth and over the point where a ray meets the structure. We test both simplifications. Direct tracing in $(\theta,\phi)$ agrees with Eq. (2) to within 0.6 percentage points (pp) (Sec. 4.3), and separate s and p channels change the conclusions by at most 0.4 pp (Table S7).

### 2.2 Reciprocity and the angular mixing number

A passive, reciprocal structure obeys detailed balance in étendue [@clarke1985]:

$$
w_i B_R(i\to j) = w_j B_R(j\to i). \tag{3}
$$

Equation (3) is the discrete form of the self-adjointness of reciprocal scattering operators [@veach1997]. Redistribution matrices of this kind are also used for textured optical sheets [@eisenlohr2015]. As a result, $M$ is similar to a symmetric matrix and its eigenvalues are real [@levin2017;@soler2026]. On the full angular range its stationary distribution is proportional to $w(1-p)$. The eigenvalue with the second-largest magnitude, $\lambda_2$, sets how fast any angular distribution relaxes to the stationary one. We express this rate as a number of encounters,

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|} , \tag{4}
$$

and call $N_{\mathrm{mix}}$ the angular mixing number. A planar interface has $N_{\mathrm{mix}} = \infty$, and a Lambertian scatterer has $N_{\mathrm{mix}} = 1$. Figures 1(b)–1(d) show $M$ for a planar interface, a hemispherical MLA, and a scattering layer. In all structures studied here $\lambda_2$ is real and positive.

Bins near grazing incidence carry almost no étendue, but they can keep their angle for many encounters and would then dominate $\lambda_2$. We therefore compute the spectral quantities ($\lambda_2$, $N_{\mathrm{mix}}$, and $p_\infty$ in Sec. 2.6) on the bins below 85°, which leave out 0.76% of the étendue. This truncated block of $B_R$ loses light to the omitted bins, so its columns sum to $s_i \le 1 - p_i$, where $s_i$ is the fraction returned into the kept bins. For $\lambda_2$ we divide each column of the block by $s_i$, which gives a column-stochastic $M$ with stationary distribution proportional to $w s$ (Note S1). For $p_\infty$ we use the block without this renormalization. Everything else, including $\eta_{\mathrm{sta}}$, $p$, $\sigma_p^2$, and $p_{\mathrm{mix}}$, is computed on the full 0–90° range. Moving the cutoff between 80° and 90° changes the mean error of the reduced model by less than 0.3 pp (Note S1).

**Fig. 1.** Multi-pass extraction as a Markov chain on angle. (a) A passive structure transmits a fraction $p(\theta)$ of the incident light and returns the rest with $B_R(\theta_i\to\theta_r)$; the device reflects the returned light with $R(\theta)$. (b)–(d) Return matrices $M$ ($\log_{10}$ scale, floor $10^{-3}$) of a planar interface, the hemispherical MLA, and a scattering layer ($g$ = 0.8, $\tau^*$ = 1). The off-diagonal bands in (c) are fixed paths of rays that are totally reflected at the lens flank, for example from 70° to 34°; they persist in a 6 × 6 supercell and with twice as many rays (Note S10). (e) Angular distribution of the power that meets the hemispherical MLA at encounter $j$, $\mathbf{v}^{(j)}$ of Eq. (2), normalized to unit area (Al device, 550 nm). (f) Escape probability at encounter $j$, $p_j$, divided by $p$, for five structures. Open symbols on the right are the complete-mixing reference values $p_{\mathrm{mix}}/p$ of Eq. (6). They are not predicted end points for slowly mixing structures such as the hemispherical MLA. The dotted line marks the fully randomized model. *(`fig1_mixing_concept.png`)*

### 2.3 Reciprocal reduced model

The reduced model keeps the $p(\theta)$ of a structure and replaces $M$ with a family that has one parameter. A fraction $1-m$ of the returned light keeps its angle, and a fraction $m$ is spread over a fixed distribution $q$:

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right]. \tag{5}
$$

Equation (3) holds only for $q \propto w(1-p)$. In words, the light spread into an angle is proportional to how much light that angle returns. The same rank-one form appears in reciprocal reflectance models that combine a diffuse lobe with an angle-dependent specular term [@ashikhmin2000;@kelemen2001]. Multi-pass extraction models, in contrast, assume isotropic internal radiance or spread the returned light isotropically [@yablonovitch1982;@schnitzer1993;@green2002;@neyts2006]. We are not aware of one that uses the reciprocal form. Isotropic spreading puts part of the returned light inside the escape cone, although the escape cone itself returns almost nothing. This breaks Eq. (3) whenever $p$ depends on angle. The kernel of Eq. (5), normalized, has the eigenvalues 1 and $1-m$, so $m = 1-\lambda_2$.

### 2.4 Complete mixing at the escape interface

For $m = 1$ the return matrix has rank one, $B_R(i\to r) = (1-p_i)\,q_r$. The returned light then escapes at the next encounter with probability

$$
p_{\mathrm{mix}} = \sum_r q_r p_r = p - \frac{\sigma_p^2}{1-p}, \qquad \sigma_p^2 = \frac{\sum_k w_k (p_k - p)^2}{\sum_k w_k}, \tag{6}
$$

with $p$ and $\sigma_p^2$ taken over the full 0–90° range. The multi-pass series then has the closed form

$$
\eta_{\mathrm{mix}} = \mathbf{p}\cdot\mathbf{v}^{(0)} + \frac{(\mathbf{p}\cdot\mathbf{u})\,[(\mathbf{1}-\mathbf{p})\cdot\mathbf{v}^{(0)}]}{1-(\mathbf{1}-\mathbf{p})\cdot\mathbf{u}}, \qquad u_k = R_k\,q_k . \tag{7}
$$

For an angle-independent reflectance $1-A'$, Eq. (7) becomes $\eta_{\mathrm{mix}} = p_0 + (1-p_0)\,p_{\mathrm{mix}}(1-A')/[p_{\mathrm{mix}} + (1-p_{\mathrm{mix}})A']$, where $p_0$ is the first-pass escape probability. This equals the fully randomized model only if $p(\theta)$ is constant.

The variance term is large for refractive surfaces. They transmit inside the critical cone and reflect everything outside it. A planar interface has $p$ = 0.40 but $p_{\mathrm{mix}}$ = 0.06. The fully randomized model gives mixed light the probability $p$, and so it overestimates the gain from mixing by an amount set by $\sigma_p^2$.

Equation (7) is a reference value for rank-one reciprocal mixing, not an upper bound. Note S1 gives a two-state example in which a reciprocal kernel of higher rank beats $\eta_{\mathrm{mix}}$ by 5 pp. Because $p_{\mathrm{mix}}$ applies at every later pass, the late-pass escape ratio of Sec. 2.6 equals $\xi = 1 - \sigma_p^2/[p(1-p)]$ for complete mixing, which is below 1 unless $p(\theta)$ is constant.

### 2.5 Mixing in transmission and at the device

The penalty in Eq. (6) arises because the escape interface itself does the mixing. Reciprocity weights redistribution in reflection by the return probability of the reflecting element, which here is $1-p$. Mixing elsewhere is weighted differently.

First, consider a lossless layer with no reflection that mixes light in transmission, an idealized volume randomizer between the device and the escape interface. Every pass through it spreads the light in proportion to $w$, whatever $p$ and $R$ are. The mixed light escapes with probability $p$, and complete mixing of this kind reproduces the fully randomized model exactly. In a lossless, index-matched substrate the position of the layer does not matter (Note S11). Real scattering layers also reflect light, which lowers their own escape probability. Section 4.1 examines how close they come to the ideal.

Second, the device can mix light in reflection, for example at a diffusely reflecting electrode. The structure keeps its matrix $B_R$, and the device reflection becomes a kernel $D$. Its columns must sum to the reflectance $R_i$, and it must satisfy $w_i D(i\to j) = w_j D(j\to i)$. The rank-one kernel that meets both conditions is

$$
D(i\to j) = \frac{R_i\,R_j\,w_j}{\sum_k w_k R_k} . \tag{8}
$$

Equation (8) is an idealized reference for complete mixing in reflection at a reciprocal device with a given $R(\theta)$. It is not a model of a specific electrode. Light mixed in this way is spread in proportion to $wR$ and escapes at the next encounter with probability

$$
p_D = \frac{\sum_k w_k R_k\,p_k}{\sum_k w_k R_k} = p + \frac{\mathrm{Cov}_w(R,p)}{\langle R\rangle_w}, \tag{9}
$$

where the covariance and the mean are étendue weighted. No variance term appears. If $R$ does not depend on angle, $p_D = p$. If the device absorbs most strongly inside the escape cone, as strong microcavities do, the covariance is negative and mixing at the device is penalized too (Note S11). A more realistic electrode that reflects part of the light specularly and part with Eq. (8) is also reciprocal; its gain grows almost linearly with the diffuse part and shows the same behavior (Table S30). Table 1 collects the three cases.

**Table 1.** Completely mixed light for three placements of the mixing. The third column is a per-encounter probability: the chance that mixed light escapes the next time it meets the escape interface. The last column is the multi-pass efficiency $\eta_{\mathrm{sta}}$ of a planar interface with this mixing (Al device, 550 nm; 0.384 without mixing, 0.798 in the fully randomized model). All quantities use the full 0–90° range.

| Where the light is mixed | Light is spread in proportion to | Escape probability per encounter of the mixed light | Equals the fully randomized model if | Multi-pass $\eta_{\mathrm{sta}}$, planar interface |
|---|---|---|---|---|
| Escape interface, in reflection | $w(1-p)$ | $p - \sigma_p^2/(1-p)$ | $p(\theta)$ is constant | 0.488 |
| Lossless volume, in transmission | $w$ | $p$ | always | 0.798 |
| Device, in reflection | $wR$ | $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$ | $R$ is uncorrelated with $p$ | 0.770 |

### 2.6 Late-pass escape ratio

The mixing number tells how fast the angular distribution relaxes but not what it relaxes to. Light that has circulated many times approaches the leading eigenvector of $B_R$, the quasi-stationary distribution of the absorbing chain [@darroch1965]. It escapes with probability $p_\infty = 1-\rho(B_R)$ per encounter, where $\rho$ is the spectral radius. This is the escape rate of an open system [@altmann2013;@ryu2006;@clauss2019]. Because $B_R$ is self-adjoint, its Rayleigh quotient with the uniform-radiance vector and its column sums [@horn2012] give

$$
\min_\theta p(\theta) \le p_\infty \le \tilde p , \tag{10}
$$

where $\tilde p$ is the probability that light with uniform radiance leaves the truncated chain, either by escape or by return beyond the cutoff ($\tilde p - p < 0.008$). We define the late-pass escape ratio as $\xi = p_\infty/\tilde p$. It equals 1 when long-lived light escapes as easily as Lambertian light, and it tends to 0 when some angles neither escape nor mix. All structures satisfy Eq. (10) numerically (Note S1). $p_\infty$ describes the structure alone, for lossless return by the device. In an OLED the surviving light is governed by $\mathrm{diag}\,\mathbf{R}\,B_R$, which also contains the parasitic loss, but this operator ranks the structures in the same order (Spearman coefficient 0.96).

---

## 3. Methods

### 3.1 Device optics

Dipole emission and stack reflectance are computed with the Chance–Prock–Silbey model [@chance1978], extended to birefringent layers [@wasey2000;@moon2015], and with transfer matrices. We use measured optical constants and treat the substrate as incoherent [@mladenovski2011]. The reference stack is Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / glass ($n$ = 1.51), with a horizontal dipole ratio of 0.865. A second stack replaces Al with Ag, and a low-loss stack follows Ref. [@kim2026]. At 550 nm, $A'$ is 0.167, 0.044, and 0.026 for the three stacks (Table S2, Fig. S1). Device optics are computed at each wavelength with dispersive optical constants. The structure matrices do not depend on wavelength, because they come from ray optics with a fixed index of 1.51. Spectrally weighted results cover 400–700 nm in 10 nm steps, weighted with the emitter spectrum times the substrate-coupled fraction $\eta_{\mathrm{sub}}(\lambda)$ (Note S10).

### 3.2 Structures and ray tracing

All structures are index matched to the substrate. The library has 38 structures, 37 of them non-planar. They are a planar interface; six non-absorbing Henyey–Greenstein [@henyey1941] scattering slabs with asymmetry parameter $g$ = 0–0.9 and reduced optical thickness $\tau^*$ = 0.25–4; nine hexagonal MLAs, which are the hemispherical and freeform optima for the total EQE and for four angular bands; 17 random MLAs on a 6 × 6 supercell with independent lenslet profiles (spread $s$), heights, or positions and radii; and five Gaussian-correlated wavy surfaces with rms slopes of 0.2–1.5, similar to the buckled surfaces used on OLEDs [@koo2010]. Ray optics is adequate for surface features much larger than the wavelength [@tang1997]; the lenslets here have radii of 10 µm.

MLAs and textures are traced with a height-field Monte Carlo ray tracer. It reproduces the LightTools return matrix of the hemispherical MLA (escape probability 0.4275 against 0.4274). The scattering slabs are traced by Monte Carlo on 0.5° bins and combined with the planar exit face, including Fresnel reflection, by the adding method [@vandehulst1980]. For two slabs, direct Monte Carlo tracing of slab and face together agrees with this to within 0.0005 in $\eta_{\mathrm{sta}}$ (Table S21). The traced matrices are used as they are, without forcing them to be symmetric. Their residual asymmetry, $\sum|F-F^{\mathsf T}|/\sum|F|$ with $F_{ji} = w_i B_R(i\to j)$, is 0.05–0.13 on 1° bins and falls to 0.002–0.007 on 15° bins, as expected for Monte Carlo noise. The planar matrix is exactly symmetric. Symmetrizing the matrices changes $N_{\mathrm{mix}}$ by at most 0.20 and $\eta_{\mathrm{sta}}$ by at most 0.04 pp (Note S1).

Numerical parameters change $N_{\mathrm{mix}}$ by at most 0.24, $\xi$ by at most 0.025, and $\eta_{\mathrm{sta}}$ by at most 0.34 pp (Tables S5 and S6, Note S10). This includes the bin width. Merging 1° bins into 2°, 3°, or 5° bins changes $p$ by less than $10^{-4}$, $p_{\mathrm{mix}}$ of the lens arrays and textures by at most 0.002, and $\eta_{\mathrm{sta}}$ by at most 0.003 up to 3° bins. For the planar interface, whose $p(\theta)$ has a sharp step at the critical angle, $p_{\mathrm{mix}}$ is 0.054, 0.056, and 0.061 for 0.5°, 1°, and 3° bins (Note S1).

The slab of Fig. 2(c) has an internal isotropic source, so its launch vector is proportional to $w$. Its planar exit face transmits with the unpolarized Fresnel coefficients, and its back reflector has an angle-independent reflectance $R$. Complete mixing at the exit face uses the rank-one kernel with $q \propto w(1-p)$ of Sec. 2.4.

### 3.3 Fabricated devices and finite-geometry model

Devices with two emitters were made on 25.4 mm glass substrates of 0.7 mm thickness. The orange devices are ITO (150 nm) / TAPC (50 nm) / TCTA (10 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / B3PyMPM (45 nm) / LiF (1 nm) / Al (100 nm). The green devices use Ir(ppy)$_2$acac in the same host, with nominal transport-layer thicknesses of 105 nm (TAPC) and 42 nm (B3PyMPM). A pixel-defining layer opens a circle of 1 mm radius. For the lens devices, a second 0.7 mm glass is attached with index-matching oil. It carries a hexagonal array of hemispherical or freeform lenslets (radius 10 µm) on a 3 × 3 mm square centered on the pixel. The substrate edges are blackened. The freeform profile differs between the two emitters (Note S13).

The EQE and the angular intensity were measured with a photodiode goniometer at 0–80° in 10° steps. The on-axis EQE is first computed as if the emission were Lambertian. The total EQE is this value times the angular factor $F = 2\int_0^{\pi/2}[I(\theta)/I(0)]\sin\theta\,d\theta$. The integral assumes that the emission is the same in all azimuths. It uses a natural cubic spline through the measured points and sets $I(90°) = 0$, because the goniometer does not reach 80–90°. The model is evaluated in exactly the same way: it computes the photodiode readings at the same angles and integrates them with the same spline, so model and measurement share any integration bias (Note S13).

The devices are modeled with the chain extended to lateral position. Between two encounters with the front face, light moves sideways by $2t\tan\theta$, where $t$ is the substrate thickness. The back-face reflectance depends on the region of the layout (pixel, Al-covered, or uncovered, each with or without ITO). Light that reaches the blackened edges is lost. The lens kernels are traced with the orientation of the hexagonal lattice resolved, and the goniometer is modeled with its scan plane, rotation axis, and detector size. For an infinite substrate this model reproduces $\eta_{\mathrm{sta}}$ to within $2\times10^{-4}$ (Note S13).

The cavity and one internal efficiency per emitter (0.89 for orange and 0.95 for green) are fixed on the bare devices. For the green stack, matching the bare angular profile required transport-layer thicknesses different from the nominal ones (65 nm B3PyMPM and 80 nm TAPC; Note S13). The lens devices are then predicted with no further parameters. The orientation of the lens lattice relative to the scan plane was not recorded, so the predictions are given as ranges over lattice orientation, scan plane, rotation axis, and photodiode size.

---

## 4. Results

### 4.1 Where the light is mixed

The angle-resolved series shows directly what repeated passes do (Al device, 550 nm). When light first reaches the hemispherical MLA, 28% of it lies beyond 60° [Fig. 1(e)]. At each pass the lenses remove light in and near the escape cone and return the rest close to its original angle. This share therefore grows to 42% after two returns and to 83% after ten, and the peak of the distribution moves from 63° to 83°. The escape probability per encounter falls with it, from 0.43 at the first encounter to 0.14 at the eleventh, a third of $p$ = 0.42 [Fig. 1(f)]. For the planar interface it drops from 0.35 to 0.05 at the second encounter. Structures that mix within one or two encounters behave differently, but they do not return to $p$ either. The scattering layer settles at 0.316 and the steep wavy surface at 0.366. Both values are within 0.003 of $p_{\mathrm{mix}}$ (0.318 and 0.369) and below $p$ (0.326 and 0.422). A random MLA, with $N_{\mathrm{mix}}$ = 1.9, lies in between. So light that a structure has mixed in reflection escapes with $p_{\mathrm{mix}}$, as Eq. (6) predicts, and not with the $p$ assumed by the fully randomized model.

The size of this penalty varies widely [Fig. 2(a)]. The ratio $p_{\mathrm{mix}}/p$ is 0.81–0.93 for lens arrays, 0.37–0.90 for wavy surfaces, and 0.14 for the planar interface. For the device stacks considered here, mixing at the device largely avoids the penalty, because their reflectance is nearly uncorrelated with $p(\theta)$ ($|p_D - p| < 0.008$, Note S11). The difference in $\eta_{\mathrm{sta}}$ between complete mixing at the device and complete mixing at the escape interface grows with $p - p_{\mathrm{mix}}$ [Fig. 2(b)]. The Pearson correlation is 0.97 over the 38 structures, 0.98 without the planar interface, and the Spearman correlation is 0.95. For the planar interface, $\eta_{\mathrm{sta}}$ is 0.384 without mixing, 0.488 with mixing at the interface, 0.770 with mixing at the device, and 0.798 with a lossless transmissive randomizer, which is the fully randomized value. For lens arrays, whose $\sigma_p^2$ is small, mixing at the device adds 1.5–4.1 pp, against at most 1.8 pp at the escape interface. Mixing at the device comes within 2.2 pp of the fully randomized model for all non-planar structures in the Al device and within 1.2 pp in the Ag and low-loss devices.

Transmissive mixing does not need an ideal randomizer. A non-absorbing Henyey–Greenstein volume under the planar exit face mixes light in transmission, but it also reflects part of it [Fig. 2(d)]. With forward scattering ($g$ = 0.8–0.9) and $\tau^*$ ≈ 0.5, it reaches $\eta_{\mathrm{sta}}$ = 0.764 in the Al device. This is 92% of the way from the planar interface (0.384) to the ideal transmissive value (0.798), and far above the 0.488 from complete reflective mixing at the interface. In the Ag and low-loss devices the same volume recovers 97% and 99%. Thicker volumes fall back to 0.58 at $\tau^*$ = 4, because they reflect more light back toward the device. The diffuse reflectance of the volume itself rises from 0.31 at $\tau^*$ = 0.5 to 0.75 at $\tau^*$ = 4, and the escape probability of the volume with its exit face drops from 0.37 to 0.19 (Table S38). Each return costs a device round trip, so efficiency falls.

The same distinction holds outside OLEDs. Figure 2(c) shows a slab with an internal isotropic source and an angle-independent back reflectance $R$ = 0.90. Here complete mixing at the back reflector is the same as complete transmissive mixing. With a planar exit face and $n$ = 1.5, the extraction efficiency is 0.44 without mixing, 0.60 with complete mixing at the exit face, and 0.87 with complete mixing at the back reflector. At $n$ = 3.5 the three values are 0.08, 0.19, and 0.37. A textured exit face transmits a similar fraction from every internal angle, so its $\sigma_p^2$ is small (Table S28). Hemispherical-MLA and steep wavy faces then come within 1.1 pp of the fully randomized result at $n \ge 2.4$, and complete mixing at the exit face changes their efficiency by at most 2.0 pp. A shallow wavy face (rms slope 0.4) stays up to 5.8 pp below the fully randomized result. Mixing at its exit face recovers none of this gap, whereas mixing at the back reflector recovers all of it.

**Fig. 2.** Where the light is mixed. (a) Escape probability of light completely mixed by the escape interface, $p_{\mathrm{mix}}$, against the escape probability $p$ of light with uniform radiance, for all structures. On the dotted line, $p_{\mathrm{mix}} = p$. (b) Difference in $\eta_{\mathrm{sta}}$ between complete mixing at the device and at the escape interface against $p - p_{\mathrm{mix}}$ (Al device, 550 nm). The legend applies to (a), (b), and (e). (c) Extraction efficiency of a slab with an internal isotropic source and an angle-independent back reflectance $R$ against its refractive index. Lines show a planar exit face without mixing, with complete mixing at the exit face, and with complete mixing at the back reflector; for constant $R$ the last case equals transmissive mixing and the fully randomized model. Symbols show textured exit faces without mixing (Table S28). The key is on the right. (d) $\eta_{\mathrm{sta}}$ of non-absorbing Henyey–Greenstein volumes under a planar exit face against the reduced optical thickness $\tau^*$ (Al device, 550 nm), compared with the ideal transmissive randomizer, complete reflective mixing at the planar interface, and the planar interface without mixing. (e) Gain in $\eta_{\mathrm{sta}}$ predicted for complete mixing by the structure at its own $p(\theta)$ (filled symbols) and by the fully randomized model (open symbols) against $N_{\mathrm{mix}}$ (Al device, 550 nm). Arrows mark fully randomized values above the axis, with their values in pp; the planar interface (10.4 and 41.4 pp) is not shown. *(`fig2_complete_mixing.png`)*

### 4.2 Gain from complete mixing by the structure

How much would each structure gain if it mixed light completely but kept its own $p(\theta)$? At 550 nm, Eq. (7) changes $\eta_{\mathrm{sta}}$ of the 37 non-planar structures by −0.4 to +1.8 pp in the Al device. The largest increase is 2.0 pp in the Ag device and 0.4 pp in the low-loss device [Fig. 2(e)]. With spectral weighting (emission peak at 572 nm), the largest increases are 1.1, 1.4, and 2.2 pp for the Al, Ag, and low-loss devices. The fully randomized model predicts more, by 0.6–2.9 pp for lens arrays, 17 pp for the shallowest wavy surface, and 31 pp for the planar interface. Mixing at the device, in comparison, adds 1.1–3.3 pp for lens arrays in the Al device with the same spectral weighting (Table S29).

The gain is larger when the cavity emits into angles that the structure transmits poorly. At 500 nm the first-pass escape probability of the hemispherical MLA is 0.36 in the Al device and 0.33 in the Ag device, well below $p$ = 0.42. Complete mixing would then raise $\eta_{\mathrm{sta}}$ by 5.7–6.4 pp (Al) and 8.4–9.4 pp (Ag), and random lenses and wavy surfaces achieve most of this increase (Table S10). Reflective mixing is therefore useful when the first-pass angular distribution does not match $p(\theta)$.

Two wavy surfaces (rms slopes 0.4 and 0.7) exceed $\eta_{\mathrm{mix}}$, by 0.42 and 0.01 pp. The second value is below the numerical uncertainty of 0.34 pp. These are the only structures in which the escape probabilities at successive encounters are anticorrelated (Note S10).

### 4.3 Accuracy and limits of the reduced model

We compared four models with the full series of Eq. (2) in 342 cases: 38 structures, three devices, and three wavelengths. The fully randomized model overestimates $\eta_{\mathrm{sta}}$ of the planar interface by up to 60 pp (Ag device, 500 nm) and that of the hemispherical MLA by 2–3 pp [Fig. S2(a)]. The opposite assumption, that every ray keeps its angle, underestimates the scattering layers [Fig. S2(b)]. A reduced model with $m = 1-\lambda_2$ and isotropic spreading fits the lens arrays but overestimates shallow wavy surfaces by up to 10.5 pp [Fig. S2(c)]. These surfaces have a small $\lambda_2$, but they return trapped light mostly to other trapped angles. Reciprocal spreading captures this without another parameter. With it, the mean absolute error falls to 0.30 pp (mean signed error −0.24 pp) and the largest error to 3.0 pp [Fig. 3(a)]. The error is not concentrated in any device or wavelength [Fig. 3(b), Table S4].

These residuals are model errors, not simulation noise. Re-tracing eight representative structures with three different random seeds changes $\eta_{\mathrm{sta}}$ by at most 0.06 pp (standard deviation) and the reduced-model error by at most 0.02 pp. For example, the error of the shallow wavy surface is −2.71 ± 0.01 pp (Table S5).

Direct ray tracing in $(\theta, \phi)$ agrees with the full series to within 0.6 pp for eight structures (Note S9, Table S19). The remaining differences vanish when the azimuth is randomized at each pass, so they come from the azimuthal correlation that the chain neglects. We also tested nine structures that were not used to develop the model: spherical caps, paraboloids, Gaussian bumps, a cone–hemisphere mixture, wavy surfaces with a ring spectrum, and two cone arrays. For seven of them the mean error is 0.48 pp and the largest 3.3 pp. For the two cone arrays the error reaches 20 pp.

The two failures point to two conditions for a single rate $1-\lambda_2$ to describe the circulating light [Fig. 3(c)]. Both can be checked from $M$ and the first-pass distribution, without evaluating the series. First, the chain must be irreducible, meaning that every angle can eventually reach every other angle. A closed set of angles gives a second eigenvalue of one; we therefore flag a chain as reducible when $\lambda_2 > 0.999$. In the 45° cone array, light within 0–4° (0.5% of the étendue) is retroreflected into the same range. This closed set gives $\lambda_2 = 1$ even though all other angles mix. Second, the slowest mode must carry a significant share of the returned light. To measure this, we expand the deviation of the first returned light, $\mathbf{u} = B_R\mathbf{v}^{(0)}$, from the stationary distribution in the eigenvectors of $M$, with amplitudes $a_\ell$. The weighted relaxation rate is then $\lambda_{\mathrm{eff}} = \sum_{\ell\ge2} a_\ell^2\lambda_\ell / \sum_{\ell\ge2} a_\ell^2$ (Note S9). The 60° cones have $\lambda_2$ = 0.985 but $\lambda_{\mathrm{eff}}$ = 0.43. All 45 other structures have $\lambda_2 - \lambda_{\mathrm{eff}} \le 0.35$ and a largest error of 3.3 pp.

We found these conditions after the out-of-sample test. We therefore tested them prospectively on twelve more structures: cones with side angles of 30–70°, square and inverted pyramids, and one-dimensional prism gratings. The predictions were recorded before the series was evaluated (Note S9). The irreducibility condition held. The 45° pyramids and the 45° prisms form closed sets, were predicted to fail, and have errors of 21 and 23 pp. The gap condition did not discriminate. It flagged two structures that the model describes to within 2.5 pp, and the eight structures that passed both conditions have a mean error of 1.4 pp but a largest error of 7.5 pp (40° cones). The reduced model is therefore less accurate for faceted structures than for the curved lenses, scatterers, and wavy surfaces of the library. Its mean error is still below that of the fully randomized model (2.8 pp). Direct traces confirm the full series for the 60° cones to within 0.6 pp, so these failures belong to the one-rate reduction, not to the angle-resolved chain.

**Fig. 3.** Accuracy of the reciprocal reduced model. (a) $\eta_{\mathrm{sta}}$ from the reduced model, with $m = 1-\lambda_2$ and $q \propto w(1-p)$, against the full series for 342 cases. (b) Error in (a) for each device (LL, low loss) and wavelength. Points are spread horizontally for visibility, and bars give the mean absolute error. (c) Largest error over the nine device–wavelength cases against $\lambda_2 - \lambda_{\mathrm{eff}}$ (Al device, 550 nm) for the 38 library structures, nine out-of-sample structures, and twelve structures tested prospectively. Crosses mark reducible chains (45° cones, pyramids, and prisms). The shaded region, gap > 0.35, is the threshold identified after the out-of-sample test. *(`fig3_reduced_model.png`)*

### 4.4 Escape probability and mixing number of the structures

The structure classes occupy separate regions of the $(p, N_{\mathrm{mix}})$ plane [Fig. 4(a)]. The planar interface has $p$ = 0.40 and $N_{\mathrm{mix}} = \infty$. Scattering layers relax within 1.1–2.2 encounters but have lower escape probabilities, $p$ = 0.19–0.39. Lens arrays have the highest escape probabilities, $p$ = 0.42–0.43, with $N_{\mathrm{mix}}$ = 6–8 for the hemispherical and 3–7.5 for the freeform arrays.

Randomizing a lens array lowers $N_{\mathrm{mix}}$ without lowering $p$ [Fig. 4(b), Table S3]. Independent lenslet profiles reduce $N_{\mathrm{mix}}$ from 7.0 to 1.9. A position jitter of 30% of the pitch reduces it to 3.5, and a height spread of ±60% only to 4.1 (Table S5). Wavy surfaces with rms slopes of 0.4–1.5 reach $N_{\mathrm{mix}}$ = 1.1–1.2, as low as a scattering layer.

The late-pass escape ratio adds information that $N_{\mathrm{mix}}$ lacks [Fig. 4(c)]. For lens arrays, $\xi$ falls from 0.80 to 0.39 as $N_{\mathrm{mix}}$ rises from 2 to 8, because light that keeps its angle collects where the lenses transmit least. The shallowest wavy surface ($\xi$ = 0.35) and a thin isotropic scattering layer ($\xi$ = 0.54) lie below this trend at $N_{\mathrm{mix}} \approx 2$. They relax quickly but return trapped light mainly to other trapped angles, and they have the largest residual errors in Fig. 3(a).

The $\eta_{\mathrm{sta}}$ contours in Fig. 4(a) are drawn for one representative angular shape, the $p(\theta)$ of the hemispherical MLA scaled to each $p$. For this shape they are nearly vertical, so the efficiency is set mainly by $p$. Other shapes of $p(\theta)$ can tilt the contours, which is one reason why $\xi$ is reported separately. We therefore checked the conclusion with each structure's own $p(\theta)$. For lenses, scattering layers, and wavy surfaces, the full range of $N_{\mathrm{mix}}$ changes $\eta_{\mathrm{sta}}$ only as much as a change of $p$ by 0.02–0.05 (Table S8). Randomizing the lenslet profiles, for example, reduces $N_{\mathrm{mix}}$ from 7.0 to 1.9 but changes $\eta_{\mathrm{sta}}$ by less than 1 pp (Al device, 550 nm). The gain from faster mixing (+0.9 pp) is offset by the change of $p(\theta)$ that comes with it (−2.2 pp, Table S9).

**Fig. 4.** Escape probability and mixing number of the structures. (a) Escape probability $p$ and angular mixing number $N_{\mathrm{mix}}$ (logarithmic axis) of all non-planar structures. $p$ and $N_{\mathrm{mix}}$ are properties of the structure alone. The contours show $\eta_{\mathrm{sta}}$ of the reduced model for the Al device at 550 nm, using the $p(\theta)$ shape of the hemispherical MLA scaled to each $p$; they illustrate one representative shape, not all structures. The hemispherical and freeform MLAs are the optima of the benchmark (Table S15). (b) $N_{\mathrm{mix}}$ (logarithmic axis) against each randomization parameter, normalized to its maximum. Dotted lines mark the hemispherical MLA and a scattering layer. (c) $\xi$ against $N_{\mathrm{mix}}$ (logarithmic axis). The planar interface ($N_{\mathrm{mix}} = \infty$, $\xi$ = 0) is drawn at the right edge. *(`fig4_structure_map.png`)*

### 4.5 Fabricated devices

The lens patch of the fabricated devices is only 3 mm wide [Fig. 5(a)], on a substrate that is 1.4 mm thick with the lens glass attached. Most light returned by the lenses lands outside the patch. There the planar face traps it until it reaches the blackened edges. The lenses therefore act mainly on first-pass light, and the experiment tests how the model redistributes this light in angle. The multi-pass gains of Secs. 4.1–4.4 are not tested by this experiment.

The model predicts the angular behavior of bare, hemispherical, and freeform devices for both emitters without adjustment [Figs. 5(b)–5(d), Table 2]. The hemispherical array concentrates the emission toward the normal. Its angular factor, the EQE divided by the EQE of a Lambertian emitter with the same on-axis intensity, is predicted as 0.62–0.68 for the orange and 0.65–0.71 for the green devices. It is measured as 0.66 ± 0.04 and 0.71 ± 0.02. The freeform arrays spread the first-pass light. For the orange freeform devices both the model and the measurement show an almost flat intensity between 30° and 50° [Fig. 5(b)]. The bare devices are close to Lambertian. All device means lie inside the predicted ranges, except that the orange freeform mean is 0.03 above its range. The mean angular profiles deviate from the center of the predicted bands by at most 0.04 rms between 10° and 80° [Figs. 5(b) and 5(c)]. Values for individual devices are listed in Tables S32 and S33.

**Table 2.** Angular factor of the fabricated devices. The model range covers lattice orientation, scan plane, rotation axis, and photodiode size for the 3 × 3 mm lens patch. Measured values are the mean ± standard deviation over devices, or over scans where only one device was available.

| Emitter | Device | Model | Measured | Sample |
|---|---|---|---|---|
| Orange | Bare | 0.99–1.02 | 1.01 ± 0.05 | 5 devices |
| Orange | Hemispherical MLA | 0.62–0.68 | 0.66 ± 0.04 | 5 devices |
| Orange | Freeform MLA | 0.92–1.04 | 1.07 ± 0.07 | 4 devices |
| Green | Bare | 1.00–1.02 | 1.02 | 2 devices |
| Green | Hemispherical MLA | 0.65–0.71 | 0.71 ± 0.02 | 1 device, 6 scans |
| Green | Freeform MLA | 0.83–0.91 | 0.87 ± 0.04 | 3 devices |

The measurements also show that on-axis metrics overstate what a small lens patch achieves [Fig. 5(e)]. The hemispherical devices have 1.22 (orange) and 1.21 (green) times the on-axis EQE of the bare devices, but only 0.81 and 0.85 times their total EQE. The lossless model shows the same contrast at a higher level, 1.39–1.55 on axis against 0.94–1.01 in total. The reason is that the absolute EQE of all lens devices is 84–91% of the lossless prediction. The importance of the lateral size of an outcoupling structure, and of the area over which light is collected, has also been noted for laboratory-size samples [@kovacic2021] and for pixelated displays [@kim2025].

A single angle-independent loss in the attached lens glass accounts for this offset. We describe it as a transmission of 0.925 each time light crosses the oil-coupled interface, fitted to the total EQE of the orange hemispherical devices. This is the only parameter fitted to lens devices. It predicts 0.89, 0.87, and 0.90 of the lossless EQE for the orange freeform, green hemispherical, and green freeform devices, against 0.90, 0.84, and 0.91 measured, and it changes the angular factors by at most 0.02 (Note S13). With this loss, the predicted on-axis and total ratios are 1.19–1.33 and 0.81–0.87 for the hemispherical devices and 0.83–1.03 and 0.84 for the freeform devices. The measured values are 1.21–1.22 and 0.81–0.85, and 0.80–1.00 and 0.85 [filled boxes in Fig. 5(e)]. The tested variations of the lens shape do not explain the offset, because those that lower the predicted on-axis gain to the measured value raise the angular factor above its measured range (Table S35). The loss value itself is inferred, not measured. A direct check, such as the transmission of an attached flat glass without lenses, was not made and would be the most useful next control. For a 10 × 10 mm patch the model predicts 1.28–1.77 times the bare EQE, but this projection relies on the same inferred loss.

**Fig. 5.** Fabricated devices and finite-geometry model. (a) Substrate layout with the pixel, the 3 × 3 mm lens patch, and the 10 × 10 mm patch used only in the model. (b), (c) Normalized angular intensity of the orange and green devices. Bands show the model range over the measurement geometries for the 3 × 3 mm patch, with the band center as a line. Symbols are measurements. In (b): one bare device, the mean and range of two freeform devices, and the mean ± standard deviation of five hemispherical devices. In (c): the mean ± standard deviation of the tilt-corrected scans of the bare (seven scans, two devices), hemispherical (six scans, one device), and freeform (twelve scans, three devices) devices. (d) Angular factor: measured mean ± standard deviation (Table 2), model range for the 3 × 3 mm patch, and model value for the 10 × 10 mm patch. (e) Ratio of lens to bare devices for the Lambertian-equivalent on-axis EQE (axial) and for the total EQE, both obtained with the angular factor of Sec. 3.3. Dashed outlines show the lossless model. Filled boxes show the model with a transmission of 0.925 per crossing of the lens-glass interface, fitted to the total EQE of the orange hemispherical devices. Measured values are the mean ± standard deviation over devices. *(`fig5_experiment.png`)*

### 4.6 Robustness

Table 3 lists the conditions under which the central result was tested. In every case, reflective mixing by the escape interface adds little, while mixing in transmission or at an angle-independent reflector recovers the fully randomized model. The exception noted above is a first-pass distribution that does not match $p(\theta)$, as at 500 nm.

**Table 3.** Robustness of the central result. Gains are changes of $\eta_{\mathrm{sta}}$ in pp for complete mixing, for the 37 non-planar structures unless stated otherwise.

| Condition | Result | Source |
|---|---|---|
| Al, Ag, and low-loss devices, 550 nm | Reflective mixing by the structure: −0.4 to +1.8 (Al), at most 2.0 (Ag) and 0.4 (low loss). Mixing at the device: +1.5 to +4.1 for lens arrays (Al) | Secs. 4.1 and 4.2 |
| Spectral weighting, 400–700 nm | Reflective mixing: at most 1.1, 1.4, and 2.2 (Al, Ag, low loss). Mixing at the device: 1.1–3.3 for lens arrays (Al) | Table S29 |
| First pass not matched to $p(\theta)$, 500 nm | Reflective mixing: 5.7–6.4 (Al) and 8.4–9.4 (Ag) for the hemispherical MLA | Table S10 |
| Exit medium with $n$ = 1.5–3.5 (slab, $R$ = 0.90) | Textured faces within 1.1 of the fully randomized model at $n \ge 2.4$; mixing at the exit face at most 2.0 | Sec. 4.1, Table S28 |
| Finite substrate (25 mm, absorbing edges) | Reflective mixing: −5.2 to +0.8. Transmissive mixing: −1.0 to +22.7 (six structures, three devices) | Fig. S10 |
| Partly diffuse electrode (diffuse fraction 0–1) | Gain nearly linear in the diffuse fraction, at most 2.3 | Table S30 |
| Numerical parameters, bin width, and reciprocity residual | $\eta_{\mathrm{sta}}$ changes by at most 0.34, $N_{\mathrm{mix}}$ by at most 0.24 | Note S1, Tables S5 and S6 |
| Separate s and p channels | At most 0.4 | Table S7 |
| Direct $(\theta,\phi)$ ray tracing | Within 0.6 of the full series | Note S9 |
| Angular cutoff of the spectral quantities, 80–90° | Mean error of the reduced model changes by less than 0.3 | Note S1 |

---

## 5. Discussion

### 5.1 Reflective and transmissive mixing

Reciprocity ties angular redistribution to escape. A surface that mixes light in reflection sends it back in proportion to how much each angle reflects. The mixed light therefore collects where escape is least likely, and its escape probability falls below $p$ by the variance term of Eq. (6). The fully randomized model instead assumes redistribution in proportion to étendue. A lossless transmissive mixer provides this, and so does an angle-independent reflector at the device. This is why faster mixing by a refractive structure does not by itself raise extraction. In the structures examined here, faster mixing usually comes with a reshaped or lower $p(\theta)$. Examples are the randomized lens arrays (Table S9) and scattering layers, whose efficiency peaks at an intermediate scattering strength (Fig. 2(d); see also Ref. [@shiang2004]). A volume with little backscattering can approach the fully randomized model. A refractive structure that mixes through its own escape interface generally stays below it when its escape probability varies with angle.

### 5.2 Limits of lens-shape optimization

A hemispherical array already has an escape probability close to the étendue limit for uniform radiance, $1/n^2$. Changing the lens profile then mainly changes $N_{\mathrm{mix}}$ and the shape of $p(\theta)$, and the efficiency contours of Fig. 4(a) are flat in these directions. The non-uniform emission of the cavity would in principle allow more. For the Al stack at 550 nm, matching the étendue of the escape cone to the brightest internal angles gives a first-pass bound of 0.540, against $1/n^2$ = 0.439 (Fig. S4, Note S4). Bounds of this kind follow from étendue conservation [@muschaweck2013;@lenef2018]. The lenses still transmit almost the same power as for uniform radiance (0.43), because they fill each outgoing direction from a wide range of internal angles. In a controlled LightTools benchmark on the device platform of Refs. [@kim2025;@kim2023], band-optimized freeform lenses beat an equally optimized hemisphere by at most 7% in one angular band and by 1.3% in total EQE (Note S12). Getting closer to the bound seems to require a surface with angle-selective transmission. For the structures examined, the fraction emitted into 40–60° also varies only from 0.32 to 0.36, close to the Lambertian value of 0.337 (Fig. S3, Table S11).

### 5.3 Implications for design

Before optimizing a structure for mixing, Eq. (6) estimates what mixing can add. For the lens arrays studied, $\sigma_p^2/(1-p)$ is 0.03–0.08 and the gain is at most about 2 pp. When the first-pass distribution matches $p(\theta)$ well, larger gains generally require placing the mixing where reciprocity does not penalize it. One option is a weakly absorbing, forward-scattering volume. With $g \ge 0.8$ and $\tau^* \approx 0.5$ it recovers 92–99% of the transmissive gain over a planar interface. The other is the device itself, if its reflectance is uncorrelated with $p(\theta)$, as it is for the three stacks here. For lens arrays in the Al device this adds 1.5–4.1 pp.

When comparing structures, rank them by $p$ first. $N_{\mathrm{mix}}$ and $\xi$ matter mainly when the first-pass distribution of the cavity does not match $p(\theta)$; then mixing can add 5–9 pp. Random MLAs and wavy surfaces with an rms slope of one or more combine a high $p$ with a low $N_{\mathrm{mix}}$ and a high $\xi$. Their spectrally weighted advantage over periodic MLAs is still only 0.2–0.5 pp in the Al and Ag devices, similar to the numerical uncertainty, and 2 pp in the low-loss device. Irregular microlens arrays are attractive for another reason: they reduce the angular color shift of microcavity devices [@lim2006]. The $(p, N_{\mathrm{mix}})$ plane can also be used to check numerical optimizers, including learning-based ones [@kim2021;@zhao2026], because it shows whether a search moves along or across the efficiency contours.

Use the one-rate reduced model only for irreducible chains, and expect errors of several percentage points for faceted structures (Sec. 4.3). The full matrix costs little more and has neither restriction.

The fabricated devices add two practical points. In the modeled 1.4 mm substrate, a lens patch must extend several substrate thicknesses beyond the pixel to extract recycled light. The model predicts that a 3 mm hemispherical patch lowers the EQE, about 6 mm gives half the gain of a 20 mm patch, and 10 mm gives 80% of it (Table S36). Other substrate thicknesses and pixel layouts will shift these numbers. The angle-integrated EQE should also be reported, because in the first-pass regime the on-axis EQE overstates the effect of a lens patch by a factor of about 1.5.

### 5.4 Scope and limitations

The analysis applies to passive external structures that are index matched to the substrate and large enough for ray optics. It does not cover internal structures that extract waveguided modes, or the change of the cavity they cause, nor angle-selective filters, absorbing scatterers, or wavelength converters. Except in Secs. 2.5 and 4.1, the device is treated as a reflector that keeps the angle, as in earlier substrate models [@greiner2007;@shiang2004]. For a partly diffuse electrode the gain grows nearly linearly with the diffuse fraction (Table S30), but the product of the device kernel with $B_R$ need not have a real spectrum. The experiment probes the first-pass regime of small lens patches, and the lens-glass loss is inferred rather than measured (Sec. 4.5). The multi-pass gains of Secs. 4.1–4.4 are computational predictions. Light trapping in solar cells is the reciprocal problem, in which a low escape probability is the goal [@deckman1983;@green2002]; we do not analyze it here.

---

## 6. Conclusion

Optical reciprocity ties the value of angular mixing to where it happens. Light mixed in reflection by the escape interface escapes with probability $p - \sigma_p^2/(1-p)$. Light mixed in transmission through a lossless volume escapes with $p$, and light mixed at the device with $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$. The fully randomized model of statistical ray optics corresponds to the second case. For 37 scattering, microlens, random-lens, and wavy structures, complete mixing by the structure raises the spectrally weighted extraction by at most 2.2 pp, while the fully randomized model predicts up to 19 pp at 550 nm. A forward-scattering volume recovers 92–99% of the transmissive gain over a planar interface. Reflective mixing helps more when the cavity emits mainly into angles that the structure transmits poorly, by 5–9 pp for the hemispherical MLA at 500 nm. These results were tested for three device stacks and for exit media with refractive indices up to 3.5, and selected finite-substrate cases are summarized in Table 3. The reciprocal reduced model predicts multi-pass extraction with a mean absolute error of 0.30 pp for curved and random structures. It fails, as predicted, for reducible chains and is less accurate for faceted structures. On fabricated devices with two emitters, the finite-geometry model reproduces the angular factors without adjustment, and it shows that small lens patches raise the on-axis EQE while lowering the total extraction. In the configurations examined, the routes predicted to raise extraction substantially are mixing in transmission or at the device, lens patches much wider than the substrate thickness, and angle-selective structures.

---

## Funding

[TBD]

## Acknowledgment

[TBD]

## Disclosures

The authors declare no conflicts of interest.

## Data availability

The per-structure arrays $p(\theta)$ and return matrices $B_R$ of the 38 library structures, the device inputs ($\mathbf{v}^{(0)}$ and $R(\theta)$ for three stacks at 400–700 nm), the source data of all figures, and Python scripts that reproduce Figs. 2–4 are not publicly available at this time but may be obtained from the authors upon reasonable request. Supplement 1 lists the scripts and data files for every figure and table.

## Supplemental document

See Supplement 1 for supporting content.
