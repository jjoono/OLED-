# Reflective and transmissive angular mixing affect substrate light extraction differently in organic light-emitting diodes

**Authors** — [TBD]

---

## Abstract

Angular mixing is widely assumed to help extract the light trapped in the substrate of organic light-emitting diodes, as in the fully randomized model of statistical ray optics. Using optical reciprocity, we show that the benefit depends on where the light is mixed. An escape interface that mixes light in reflection returns it mainly to angles from which it cannot escape, so the mixed light escapes less often, by an amount set by the angular variance of the escape probability. Mixing in transmission through a lossless volume has no such penalty and reproduces the fully randomized model. For 37 microlens, scattering, and textured structures, complete mixing by the structure adds at most 2.2 percentage points to the extraction efficiency, whereas a forward-scattering volume recovers 92–99% of the fully randomized gain. Microlens devices with two emitters confirm the predicted angular redistribution.

---

## 1. Introduction

Organic light-emitting diodes (OLEDs) are widely used in displays and lighting. Their internal quantum efficiency is close to unity, so their efficiency is now limited mainly by light extraction [@brutting2013]. A planar OLED emits only part of the light it generates. The rest is lost to waveguided and surface-plasmon modes, to electrode absorption, and to total internal reflection at the substrate–air interface [@brutting2013;@greiner2007]. To recover the light trapped in the substrate, researchers have added microlens arrays (MLAs) [@moller2002;@sun2006;@qu2018], scattering layers [@song2018;@levell2013], or textured films [@koo2010;@kovacic2022] to the substrate, and such structures are now also designed by machine learning [@kim2021;@zhao2026]. Even so, the reported external quantum efficiencies (EQEs) mostly stay between 50% and 60% [@song2018;@kim2025;@kim2023;@li2019], because unextracted light is partly absorbed each time it returns to the device [@kim2026].

Two kinds of model describe this multi-pass process. Statistical ray optics assumes that the light is fully randomized at every pass [@yablonovitch1982;@schnitzer1993;@green2002]. The substrate-to-air efficiency is then $\eta = p/[p+(1-p)A']$, where $p$ is the escape probability of light with uniform radiance and $A'$ is the round-trip loss [@kim2026]. Angle-resolved models follow the angular distribution from pass to pass with redistribution matrices [@greiner2007;@shiang2004;@neyts2006;@eisenlohr2015;@david2013;@altazin2015;@kovacic2021]. They are accurate but give no rule for comparing structures. The two descriptions differ most for the structures of practical interest. A hemispherical MLA returns most unextracted light close to its angle of arrival, so the surviving light collects at angles that the lenses transmit poorly, and the fully randomized model overestimates extraction [@kim2026]. This suggests designing structures that mix angles faster, as studied for textured absorbers [@deckman1983;@manzoor2020;@battaglia2012]. However, whether faster mixing actually raises extraction has not been shown.

In this work, we show that the benefit of angular mixing depends on where the light is mixed, and that optical reciprocity fixes this dependence. An escape interface that mixes light in reflection must return it to each angle in proportion to how much that angle reflects. Completely mixed light then escapes with probability $p - \sigma_p^2/(1-p)$, where $\sigma_p^2$ is the angular variance of the escape probability. Mixing in transmission through a lossless volume has no such penalty, and mixing at the device is penalized only if the device reflectance is correlated with the escape probability. The fully randomized model is therefore a model of transmissive mixing and does not describe the reflective mixing that refractive outcoupling structures provide. We also derive a reduced model with one mixing parameter and evaluate 38 structures in three device stacks. Complete mixing by the structure itself adds at most 2.2 percentage points (pp) to the extraction efficiency, whereas a forward-scattering volume recovers most of the fully randomized gain. Finally, we compare a finite-geometry version of the model with fabricated MLA devices of two emitters.

---

## 2. Theory and methods

### 2.1 Angle-resolved model and reciprocity

Figure 1(a) illustrates the model. We divide the polar angle inside the substrate into 1° bins $k$ with étendue $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$, so that uniform radiance means power proportional to $w_k$. Each time light meets the structure, a fraction $p_k$ of the light in bin $k$ leaves into air, and a fraction $B_R(i\to r)$ of the light in bin $i$ returns into bin $r$. The structures do not absorb, so $p_i + \sum_r B_R(i\to r) = 1$. We write the return matrix as

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}), \tag{1}
$$

where each column of $M$ sums to one. $M$ is thus the transition matrix of a Markov chain on angle, and $p = \sum_k w_k p_k/\sum_k w_k$ is the escape probability of light with uniform radiance.

The launch vector $\mathbf{v}^{(0)}$ holds the power that the emitter couples into each bin, summed over azimuth and polarization and normalized to $\sum_k v^{(0)}_k = 1$. The device reflects returned light with its angle-resolved reflectance $R(\theta)$ and, unless stated otherwise, keeps its angle. If $\mathbf{v}^{(j)}$ is the light that meets the structure for the $(j+1)$th time, then

$$
\mathbf{v}^{(j+1)} = \mathrm{diag}\,\mathbf{R}\; B_R\,\mathbf{v}^{(j)} , \qquad \eta_{\mathrm{sta}} = \sum_{j\ge 0}\mathbf{p}^{\mathsf T}\mathbf{v}^{(j)} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} , \tag{2}
$$

where $\eta_{\mathrm{sta}}$ is the substrate-to-air extraction efficiency [@kim2026]. The escape probability at encounter $j$ is $p_j = \mathbf{p}^{\mathsf T}\mathbf{v}^{(j)}/\sum_k v^{(j)}_k$. The scalar $A' = 1-\sum_k w_k R_k/\sum_k w_k$ is used only for comparison with the fully randomized model.

A passive, reciprocal structure obeys detailed balance in étendue [@clarke1985;@veach1997]:

$$
w_i B_R(i\to j) = w_j B_R(j\to i). \tag{3}
$$

$M$ is therefore similar to a symmetric matrix, and its eigenvalues are real [@horn2012]. The eigenvalue with the second-largest magnitude, $\lambda_2$, sets how fast an angular distribution relaxes to the stationary one. We express this rate as a number of encounters,

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|} , \tag{4}
$$

and call $N_{\mathrm{mix}}$ the angular mixing number. It is the per-pass analog of the transport factor in radiative transfer [@vanrossum1999] and of the relaxation time of a reversible Markov chain [@levin2017], and similar spectral descriptions are established for light-transport operators and open systems [@soler2026;@koenderink1983;@darroch1965;@altmann2013;@ryu2006;@clauss2019]. Figures 1(b)–1(d) show $M$ for three structures. A planar interface keeps every angle and has $N_{\mathrm{mix}} = \infty$, whereas a Lambertian scatterer has $N_{\mathrm{mix}} = 1$. The largest eigenvalue of $B_R$ describes light that has circulated many times. Such light escapes with probability $p_\infty = 1-\rho(B_R)$ per encounter, where $\rho$ is the spectral radius. We define the late-pass escape ratio $\xi = p_\infty/\tilde p$, where $\tilde p$ is the corresponding probability for light with uniform radiance ($\tilde p - p < 0.008$). It equals 1 when long-lived light escapes as easily as Lambertian light and tends to 0 when some angles neither escape nor mix. The spectral quantities are computed on the bins below 85° (Supplement 1, Note S1).

**Fig. 1.** Multi-pass extraction as a Markov chain on angle. (a) The structure transmits a fraction $p(\theta)$ of the incident light and returns the rest with $B_R$; the device reflects it with $R(\theta)$. (b)–(d) Return matrices $M$ ($\log_{10}$ scale) of (b) a planar interface, (c) the hemispherical MLA, and (d) a scattering layer ($g$ = 0.8, $\tau^*$ = 1). (e) Angular distribution of the light that meets the hemispherical MLA at encounter $j$ (Al device, 550 nm). (f) Escape probability at encounter $j$ relative to $p$. Open symbols: $p_{\mathrm{mix}}/p$ from Eq. (6). Dotted line: fully randomized model. *(`fig1_mixing_concept.png`)*

### 2.2 Reduced model and complete mixing at the escape interface

The reduced model keeps the $p(\theta)$ of a structure and replaces $M$ with a one-parameter family. A fraction $1-m$ of the returned light keeps its angle, and a fraction $m$ is spread over a fixed distribution $q$:

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right]. \tag{5}
$$

Equation (3) holds only for $q \propto w(1-p)$: the light spread into an angle is proportional to how much light that angle returns. The same form appears in reciprocal reflectance models [@ashikhmin2000;@kelemen2001]. Earlier multi-pass extraction models instead spread the returned light isotropically [@yablonovitch1982;@schnitzer1993;@neyts2006], which violates Eq. (3) whenever $p$ depends on angle. The kernel of Eq. (5) has the eigenvalues 1 and $1-m$, so we set $m = 1-\lambda_2$.

For $m = 1$ the structure mixes the returned light completely. This light escapes at the next encounter with probability

$$
p_{\mathrm{mix}} = \sum_r q_r p_r = p - \frac{\sigma_p^2}{1-p}, \qquad \sigma_p^2 = \frac{\sum_k w_k (p_k - p)^2}{\sum_k w_k}, \tag{6}
$$

and the multi-pass series has the closed form

$$
\eta_{\mathrm{mix}} = \mathbf{p}\cdot\mathbf{v}^{(0)} + \frac{(\mathbf{p}\cdot\mathbf{u})\,[(\mathbf{1}-\mathbf{p})\cdot\mathbf{v}^{(0)}]}{1-(\mathbf{1}-\mathbf{p})\cdot\mathbf{u}}, \qquad u_k = R_k\,q_k . \tag{7}
$$

Equation (7) equals the fully randomized model only if $p(\theta)$ is constant. The variance term is large for refractive surfaces, which transmit inside the critical cone and reflect everything outside it. A planar interface has $p$ = 0.40 but $p_{\mathrm{mix}}$ = 0.06. The fully randomized model gives mixed light the probability $p$ and so overestimates the gain from mixing. Equation (7) is a reference value for rank-one reciprocal mixing, not an upper bound (Note S1).

### 2.3 Mixing in transmission and at the device

The penalty in Eq. (6) arises because the escape interface itself does the mixing. Reciprocity weights redistribution in reflection by the return probability of the reflecting element, which here is $1-p$. Mixing elsewhere is weighted differently (Table 1).

A lossless, non-reflecting layer that mixes light in transmission spreads it in proportion to $w$ at every pass. The mixed light escapes with probability $p$, and the fully randomized model is recovered exactly. In an index-matched substrate the position of the layer does not matter (Note S11).

The device can also mix light in reflection, for example at a diffusely reflecting electrode. Its reflection then becomes a kernel $D$ whose columns sum to $R_i$ and which satisfies $w_i D(i\to j) = w_j D(j\to i)$. The rank-one kernel that meets both conditions is

$$
D(i\to j) = \frac{R_i\,R_j\,w_j}{\sum_k w_k R_k} . \tag{8}
$$

It spreads the light in proportion to $wR$, and the mixed light escapes at the next encounter with probability

$$
p_D = \frac{\sum_k w_k R_k\,p_k}{\sum_k w_k R_k} = p + \frac{\mathrm{Cov}_w(R,p)}{\langle R\rangle_w}, \tag{9}
$$

where the covariance and the mean are étendue weighted. No variance term appears. If $R$ does not depend on angle, $p_D = p$. If the device absorbs most strongly inside the escape cone, as strong microcavities do, the covariance is negative and mixing at the device is penalized as well (Note S11).

**Table 1.** Escape probability of completely mixed light for three placements of the mixing. The last column is $\eta_{\mathrm{sta}}$ of a planar interface in the Al device at 550 nm (0.384 without mixing).

| Where the light is mixed | Light is spread in proportion to | Escape probability of the mixed light | $\eta_{\mathrm{sta}}$, planar interface |
|---|---|---|---|
| Escape interface, in reflection | $w(1-p)$ | $p - \sigma_p^2/(1-p)$ | 0.488 |
| Lossless volume, in transmission | $w$ | $p$ | 0.798 |
| Device, in reflection | $wR$ | $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$ | 0.770 |

### 2.4 Optical simulation

Dipole emission and stack reflectance are computed with the Chance–Prock–Silbey model [@chance1978], extended to birefringent layers [@wasey2000;@moon2015], and with transfer matrices. We use measured optical constants and treat the substrate as incoherent [@mladenovski2011]. The reference stack is Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / glass ($n$ = 1.51). A second stack replaces Al with Ag, and a low-loss stack follows Ref. [@kim2026]. Their $A'$ at 550 nm is 0.167, 0.044, and 0.026 (Table S2). Results are given at 550 nm or spectrally weighted over 400–700 nm with the emitter spectrum (Note S10).

The library has 38 structures, all index matched to the substrate: a planar interface; six non-absorbing Henyey–Greenstein [@henyey1941] scattering layers with asymmetry parameter $g$ = 0–0.9 and reduced optical thickness $\tau^*$ = 0.25–4; nine hexagonal MLAs with hemispherical or freeform lenslets; 17 random MLAs with independent lenslet profiles, heights, or positions; and five Gaussian wavy surfaces with rms slopes of 0.2–1.5 (Note S2). Lens arrays and textures are traced with a height-field Monte Carlo ray tracer, which reproduces the LightTools escape probability of the hemispherical MLA to within $10^{-4}$. Ray optics is adequate for the 10 µm lenslets [@tang1997]. The scattering layers are traced by Monte Carlo and combined with the planar exit face by the adding method [@vandehulst1980]. The traced matrices are used without symmetrization. Bin width, ray number, the residual asymmetry of the matrices, and polarization change $\eta_{\mathrm{sta}}$ by at most 0.4 pp (Note S1, Tables S5–S7).

### 2.5 Device fabrication and characterization

Orange and green bottom-emitting devices were made on 25.4 mm glass substrates of 0.7 mm thickness, with a circular pixel of 1 mm radius. The orange stack is ITO (150 nm) / TAPC (50 nm) / TCTA (10 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / B3PyMPM (45 nm) / LiF (1 nm) / Al (100 nm). The green devices use Ir(ppy)$_2$acac in the same host (Note S13). For the lens devices, a second 0.7 mm glass is attached with index-matching oil. It carries a hexagonal array of hemispherical or freeform lenslets (radius 10 µm) on a 3 × 3 mm patch centered on the pixel. The substrate edges are blackened.

The angular intensity $I(\theta)$ was measured with a photodiode goniometer from 0° to 80° in 10° steps. The total EQE is the on-axis EQE of a Lambertian emitter times the angular factor $F = 2\int_0^{\pi/2}[I(\theta)/I(0)]\sin\theta\,d\theta$, integrated with a spline through the measured points.

The devices are modeled with the chain extended to lateral position. Between two encounters with the front face, light moves sideways by $2t\tan\theta$, where $t$ is the substrate thickness. The back reflectance follows the device layout, and light that reaches the edges is lost (Note S13). The cavity and one internal efficiency per emitter are calibrated on the bare devices. The lens devices are then predicted without further parameters, and their angular factors are integrated in the same way as the measured ones. Model values are ranges over lens-lattice orientation and goniometer geometry.

---

## 3. Results and discussion

### 3.1 Reflective and transmissive mixing

Figure 1(e) shows the angular distribution of the light that meets the hemispherical MLA at successive encounters (Al device, 550 nm). At the first encounter, 28% of the light lies beyond 60°. After ten returns this share is 83%, and the escape probability per encounter has fallen to a third of $p$ [Fig. 1(f)]. Repeated passes thus concentrate the light where it cannot escape. Structures that mix within one or two encounters do not return to $p$ either. A scattering layer and a steep wavy surface settled within 0.003 of their $p_{\mathrm{mix}}$, as predicted by Eq. (6).

Figure 2(a) compares $p_{\mathrm{mix}}$ with $p$ for all structures. The ratio $p_{\mathrm{mix}}/p$ is 0.81–0.93 for lens arrays, 0.37–0.90 for wavy surfaces, and 0.14 for the planar interface. Mixing at the device avoids this penalty in the three stacks studied, because their reflectance is nearly uncorrelated with $p(\theta)$ ($|p_D - p| < 0.008$). As shown in Fig. 2(b), the advantage of mixing at the device over mixing at the escape interface therefore grows with $p - p_{\mathrm{mix}}$. It is largest for the planar interface (Table 1) and amounts to a few pp for lens arrays.

Figure 2(e) shows the gain that each structure would obtain if it mixed the light completely but kept its own $p(\theta)$. For the 37 non-planar structures, Eq. (7) changed $\eta_{\mathrm{sta}}$ by −0.4 to +1.8 pp in the Al device at 550 nm, and by at most 2.2 pp with spectral weighting in any of the three devices. The fully randomized model predicts up to 19 pp for the same structures. Reflective mixing helps more when the cavity emits mainly into angles that the structure transmits poorly. At 500 nm, where the first-pass escape probability of the hemispherical MLA is well below $p$, complete mixing would add 6–9 pp, and random lenses and wavy surfaces realize most of this gain (Table S10).

Figure 2(d) shows that transmissive mixing does not require an ideal randomizer. A forward-scattering volume ($g$ = 0.8–0.9, $\tau^*$ ≈ 0.5) under the planar exit face reached $\eta_{\mathrm{sta}}$ = 0.764 in the Al device. This is 92% of the way from the planar interface to the ideal transmissive value, and the same volume recovered 97% and 99% in the Ag and low-loss devices. Thicker volumes lose this advantage because they reflect more light back to the device (Table S38).

The distinction is not specific to OLEDs. Figure 2(c) shows the extraction efficiency of a slab with an internal isotropic source and a back reflectance of 0.90. At $n$ = 1.5, complete mixing at a planar exit face raised the efficiency from 0.44 to 0.60, whereas mixing at the back reflector, which here equals transmissive mixing, gave 0.87. The ordering persists up to $n$ = 3.5. Textured exit faces have a small $\sigma_p^2$ and come close to the fully randomized result without any mixing (Table S28). The same ordering appears in light trapping, the reciprocal problem [@deckman1983;@green2002]. In a weakly absorbing slab, mixing only the light that the front face reflects never gave a longer path than a Lambertian back reflector (Table S39).

**Fig. 2.** Where the light is mixed. (a) Escape probability $p_{\mathrm{mix}}$ of light completely mixed by the escape interface against $p$. On the dotted line, $p_{\mathrm{mix}} = p$. (b) Difference in $\eta_{\mathrm{sta}}$ between complete mixing at the device and at the escape interface against $p - p_{\mathrm{mix}}$. The legend applies to (a), (b), and (e). (c) Extraction efficiency of a slab with an internal isotropic source and back reflectance $R$ against its refractive index. Lines: planar exit face without mixing, with complete mixing at the exit face, and with complete mixing at the back reflector. Symbols: textured exit faces without mixing. (d) $\eta_{\mathrm{sta}}$ of Henyey–Greenstein volumes under a planar exit face against the reduced optical thickness $\tau^*$. (e) Gain in $\eta_{\mathrm{sta}}$ for complete mixing by the structure (filled symbols) and in the fully randomized model (open symbols) against $N_{\mathrm{mix}}$. Arrows mark values above the axis. Panels (a), (b), (d), and (e) are for the Al device at 550 nm. *(`fig2_complete_mixing.png`)*

### 3.2 Reduced model and structure map

Figure 3(a) compares the reduced model with the full series of Eq. (2) for 342 cases: 38 structures, three devices, and three wavelengths. With $m = 1-\lambda_2$ and reciprocal spreading, the mean absolute error was 0.30 pp and the largest error 3.0 pp, with no concentration in any device or wavelength [Fig. 3(b)]. For comparison, the fully randomized model erred by up to 60 pp for the planar interface, and a reduced model with isotropic spreading by up to 10.5 pp for shallow wavy surfaces (Fig. S2). Direct ray tracing in polar and azimuthal angle agreed with the full series to within 0.6 pp (Note S9).

Figure 3(c) shows the limits of the reduced model, which we tested on 21 further structures chosen after the model was fixed (Note S9). The model fails when the chain is reducible, that is, when a set of angles is closed under reflection and $\lambda_2 > 0.999$. This is the case for 45° cones, pyramids, and prisms, which retroreflect near-normal light, and the errors were about 20 pp. For the pyramids and prisms, the failure was predicted from $M$ before the series was evaluated. A second diagnostic compares $\lambda_2$ with the relaxation rate weighted by the first returned light, $\lambda_{\mathrm{eff}}$ (Note S9). It identified the 60° cones, in which the slowest mode carries little power, but was not reliable in the prospective test. For the remaining faceted structures the error reached 7.5 pp. The full matrix has none of these restrictions.

Figure 4(a) shows the structures in the $(p, N_{\mathrm{mix}})$ plane, where the structure classes occupy separate regions. Lens arrays have the highest escape probabilities ($p$ = 0.42–0.43) with $N_{\mathrm{mix}}$ = 3–8. Scattering layers relax within about two encounters but have a lower $p$ (0.19–0.39). Randomizing the lenslet profiles lowered $N_{\mathrm{mix}}$ from 7.0 to 1.9 without lowering $p$ [Fig. 4(b)], and $\xi$ rose from 0.39 to 0.80 as $N_{\mathrm{mix}}$ fell from 8 to 2 [Fig. 4(c)]. However, $\eta_{\mathrm{sta}}$ changed by less than 1 pp, because the gain from faster mixing is offset by the accompanying change of $p(\theta)$ (Table S9). The efficiency is therefore set mainly by $p$, as indicated by the nearly vertical contours in Fig. 4(a).

Figure 4(a) also explains why optimizing the lens shape gives small gains. A hemispherical array already has an escape probability close to the étendue limit for uniform radiance, $1/n^2$ = 0.439. The non-uniform emission of the cavity would allow a first-pass value of up to 0.540 (Note S4) [@muschaweck2013;@lenef2018], but refractive lenslets fill each outgoing direction from a wide range of internal angles and transmit about 0.43 in either case. In a LightTools benchmark on the device platform of Refs. [@kim2025;@kim2023], optimized freeform lenses exceeded the optimized hemisphere by 1.3% in total EQE (Note S12). Approaching the bound would require a surface with angle-selective transmission.

**Fig. 3.** Accuracy of the reciprocal reduced model. (a) $\eta_{\mathrm{sta}}$ from the reduced model against the full series for 342 cases. (b) Error in (a) for each device (LL, low loss) and wavelength. Bars give the mean absolute error. (c) Largest error against $\lambda_2 - \lambda_{\mathrm{eff}}$ for the 38 library structures and 21 additional structures. Crosses mark reducible chains. The shaded region marks $\lambda_2 - \lambda_{\mathrm{eff}} > 0.35$. *(`fig3_reduced_model.png`)*

**Fig. 4.** Escape probability and mixing number of the structures. (a) $p$ and $N_{\mathrm{mix}}$ of all non-planar structures. Contours show $\eta_{\mathrm{sta}}$ of the reduced model (Al device, 550 nm) for the $p(\theta)$ shape of the hemispherical MLA scaled to each $p$. (b) $N_{\mathrm{mix}}$ against each randomization parameter, normalized to its maximum. Dotted lines mark the hemispherical MLA and a scattering layer. (c) $\xi$ against $N_{\mathrm{mix}}$. *(`fig4_structure_map.png`)*

### 3.3 Fabricated devices

Figure 5(a) shows the layout of the fabricated devices. The lens patch is 3 mm wide, on a substrate that is 1.4 mm thick with the lens glass attached. Most light returned by the lenses lands outside the patch, so the lenses act mainly on first-pass light. The experiment therefore tests how the model redistributes this light in angle.

Figures 5(b)–5(d) and Table 2 compare the measured angular behavior with the model. The model predicted the bare, hemispherical, and freeform devices of both emitters without adjustment. The hemispherical array concentrated the emission toward the normal and lowered the angular factor to about 0.7. The freeform arrays spread the first-pass light and kept the factor near 0.9–1.1. All measured means lie inside the predicted ranges, except that the orange freeform mean is 0.03 above its range. The mean angular profiles deviate from the centers of the predicted bands by at most 0.04 rms.

**Table 2.** Angular factor of the fabricated devices. Model: range over lens-lattice orientation and goniometer geometry. Measured: mean ± standard deviation over devices, or over scans where only one device was available.

| Emitter | Device | Model | Measured | Sample |
|---|---|---|---|---|
| Orange | Bare | 0.99–1.02 | 1.01 ± 0.05 | 5 devices |
| Orange | Hemispherical MLA | 0.62–0.68 | 0.66 ± 0.04 | 5 devices |
| Orange | Freeform MLA | 0.92–1.04 | 1.07 ± 0.07 | 4 devices |
| Green | Bare | 1.00–1.02 | 1.02 | 2 devices |
| Green | Hemispherical MLA | 0.65–0.71 | 0.71 ± 0.02 | 1 device, 6 scans |
| Green | Freeform MLA | 0.83–0.91 | 0.87 ± 0.04 | 3 devices |

Figure 5(e) compares the on-axis and total EQE of the lens devices with those of the bare devices. Because the hemispherical array narrows the emission, its on-axis gain overstates the extraction. The on-axis EQE was 1.2 times that of the bare devices, but the total EQE was 0.81–0.85 times. The measured EQEs of the lens devices were 84–91% of the lossless prediction. A single angle-independent transmission of 0.925 per crossing of the oil-coupled lens glass, fitted to the orange hemispherical devices, reproduced this offset for the other three device types to within 0.03 and changed the angular factors by at most 0.02 (Note S13).

The model also shows how wide a lens patch must be. To extract recycled light, it must extend several substrate thicknesses beyond the pixel. A 10 mm patch gives 80% of the gain of a 20 mm patch, whereas a 3 mm patch lowers the total EQE (Table S36). This agrees with earlier observations on laboratory samples and pixelated displays [@kovacic2021;@kim2025]. The angle-integrated EQE should therefore be reported together with the on-axis value.

**Fig. 5.** Fabricated devices and finite-geometry model. (a) Substrate layout with the pixel and the lens patch. (b), (c) Normalized angular intensity of the (b) orange and (c) green devices. Bands show the model range and symbols the measurements (mean with standard deviation or range; Tables S32 and S33). (d) Angular factor: measurement, model range for the 3 × 3 mm patch, and model value for a 10 × 10 mm patch. (e) Ratio of lens to bare devices for the on-axis and total EQE: measurement, lossless model (dashed outlines), and model with a transmission of 0.925 per crossing of the lens-glass interface (filled boxes). *(`fig5_experiment.png`)*

### 3.4 Design implications

These results suggest the following guidelines. Before a structure is optimized for mixing, Eq. (6) estimates what mixing can add. For the lens arrays studied, $\sigma_p^2/(1-p)$ is 0.03–0.08, and the gain is about 2 pp at most. Larger gains require mixing where reciprocity does not penalize it: in a forward-scattering volume, or at a device whose reflectance is uncorrelated with $p(\theta)$. The latter adds 1.5–4.1 pp for lens arrays in the Al device. Structures should be ranked by $p$ first, because $N_{\mathrm{mix}}$ and $\xi$ matter mainly when the first-pass distribution does not match $p(\theta)$. Random MLAs and steep wavy surfaces combine a high $p$ with fast mixing, but their spectrally weighted advantage over periodic MLAs is at most 2 pp. They remain attractive because they reduce the angular color shift [@lim2006]. These conclusions also hold with spectral weighting, for exit media with refractive indices up to 3.5, for finite substrates, and for partly diffuse electrodes (Table S40).

---

## 4. Conclusion

In this work, we showed that optical reciprocity ties the value of angular mixing to where it happens. Light mixed in reflection by the escape interface escapes with probability $p - \sigma_p^2/(1-p)$, light mixed in transmission through a lossless volume with $p$, and light mixed at the device with $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$. The fully randomized model of statistical ray optics corresponds to the second case. For 37 scattering, microlens, and textured structures, complete mixing by the structure added at most 2.2 pp to the spectrally weighted extraction, whereas a forward-scattering volume recovered 92–99% of the fully randomized gain over a planar interface. A one-parameter reciprocal model reproduced the full calculation with a mean error of 0.30 pp for curved and random structures, and fabricated devices with two emitters confirmed the predicted angular redistribution of first-pass light. It should be noted that the analysis applies to passive external structures that are index matched to the substrate and large enough for ray optics. We believe that these results will help in designing outcoupling structures: substantially higher extraction calls for mixing in transmission or at the device, lens patches much wider than the substrate thickness, or angle-selective structures.

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
