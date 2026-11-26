# Reflective and transmissive angular mixing affect substrate light extraction differently in organic light-emitting diodes

**Authors** — [TBD]

---

## Abstract

Angular mixing is widely assumed to help extract the light trapped in the substrate of organic light-emitting diodes (OLEDs), as in the fully randomized (FR) model of statistical ray optics. In this work, we use optical reciprocity to show that the benefit depends on where the light is mixed. An escape interface that mixes light in reflection returns it mainly to angles from which it cannot escape, so that the mixed light escapes less often, by an amount set by the angular variance of the escape probability. In contrast, mixing in transmission through a lossless volume has no such penalty and reproduces the FR model. For 37 microlens, scattering, and textured structures in three device stacks, complete mixing by the structure itself adds at most 2.2 percentage points to the light extraction efficiency, whereas a forward-scattering volume recovers 92% to 99% of the FR gain over a planar interface. Angle-resolved measurements on microMLA devices with two emitters confirm the predicted angular redistribution, and the results provide guidelines for placing angular mixing in outcoupling structures.

---

## 1. Introduction

Organic light-emitting diodes (OLEDs) are widely used in displays and lighting, and because their internal quantum efficiency is close to unity, their efficiency is now limited mainly by the light extraction efficiency (LEE). A planar OLED emits only part of the light it generates, while the rest is lost to waveguided and surface-plasmon modes, to electrode absorption, and to total internal reflection (TIR) at the interface between the substrate and air [1,2]. To recover the light trapped in the substrate by TIR, researchers have attached microlens arrays (MLAs), scattering layers, or textured films to the substrate [3–7]. However, the external quantum efficiencies (EQEs) reported with such structures mostly remain between 50% and 60% [8,9], because the unextracted light is partly absorbed each time it returns to the OLED.

Two kinds of models have been used to describe this multi-pass process. Statistical ray optics assumes fully randomized (FR) light at every pass [10–12]. The substrate-to-air efficiency is then $\eta = p/[p+(1-p)A']$, where $p$ is the escape probability of light with uniform radiance and $A'$ is the round-trip loss [13]. Angle-resolved models instead follow the angular distribution from pass to pass with redistribution matrices [14–17]; they are accurate but give no simple rule for comparing structures. The two descriptions differ most for the structures of practical interest. A hemispherical MLA returns most of the unextracted light close to its angle of arrival, and the surviving light collects at angles that the MLA transmits poorly, so the FR model overestimates the LEE. This observation suggests designing structures that mix angles faster, as has been studied for textured absorbers [18]. However, whether faster mixing actually raises the LEE has not been shown.

In this work, we show that the benefit of angular mixing depends on where the light is mixed, and that optical reciprocity fixes this dependence. An escape interface that mixes the light in reflection must return it to each angle in proportion to the reflectance of this angle. The completely mixed light then escapes with a probability of $p - \sigma_p^2/(1-p)$, where $\sigma_p^2$ is the angular variance of the escape probability. In contrast, mixing in transmission through a lossless volume has no such penalty, and mixing at the device is penalized only if the device reflectance is correlated with the escape probability. The FR model is therefore a model of transmissive mixing and does not describe the reflective mixing that refractive outcoupling structures provide. We also derive a reduced model with a single mixing parameter and evaluate 38 structures in three device stacks. Complete mixing by the structure itself adds at most 2.2 percentage points (pp) to the LEE, whereas a forward-scattering volume recovers most of the FR gain. Finally, we compare a finite-geometry version of the model with fabricated MLA devices of two emitters.

---

## 2. Theory and methods

### 2.1 Angle-resolved model and reciprocity

Figure 1(a) illustrates the model. The polar angle inside the substrate is divided into bins $k$ with a width of 1° and an étendue of $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$, so that uniform radiance corresponds to a power proportional to $w_k$. Each time the light meets the structure, a fraction $p_k$ of the light in bin $k$ leaves into air and a fraction $B_R(i\to r)$ of the light in bin $i$ returns into bin $r$. Because the structures do not absorb, $p_i + \sum_r B_R(i\to r) = 1$. We write the return matrix as

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}), \tag{1}
$$

where each column of $M$ sums to one, so that $M$ is the transition matrix of a Markov chain on angle and $p = \sum_k w_k p_k/\sum_k w_k$ is the escape probability of light with uniform radiance.

The launch vector $\mathbf{v}^{(0)}$ contains the power that the emitter couples into each bin, summed over azimuth and polarization and normalized to unity. The device reflects the returned light with its angle-resolved reflectance $R(\theta)$ and, unless stated otherwise, keeps its angle. If $\mathbf{v}^{(j)}$ denotes the light that meets the structure for the $(j+1)$th time, then

$$
\mathbf{v}^{(j+1)} = \mathrm{diag}\,\mathbf{R}\; B_R\,\mathbf{v}^{(j)} , \qquad \eta_{\mathrm{sta}} = \sum_{j\ge 0}\mathbf{p}^{\mathsf T}\mathbf{v}^{(j)} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} , \tag{2}
$$

where $\eta_{\mathrm{sta}}$ is the substrate-to-air LEE and $p_j = \mathbf{p}^{\mathsf T}\mathbf{v}^{(j)}/\sum_k v^{(j)}_k$ is the escape probability at encounter $j$. The scalar loss $A' = 1-\sum_k w_k R_k/\sum_k w_k$ is used only for comparison with the FR model.

A passive, reciprocal structure obeys detailed balance in étendue [19,20],

$$
w_i B_R(i\to j) = w_j B_R(j\to i), \tag{3}
$$

so that $M$ is similar to a symmetric matrix and has real eigenvalues. The eigenvalue with the second-largest magnitude, $\lambda_2$, determines how fast an angular distribution relaxes to the stationary one, and we express this rate as a number of encounters,

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|} , \tag{4}
$$

which we call the angular mixing number. It is the per-pass analog of the transport factor in radiative transfer and of the relaxation time of a reversible Markov chain [21,22]. Figures 1(b)–1(d) show $M$ for three structures, where a planar interface keeps every angle and has $N_{\mathrm{mix}} = \infty$, whereas a Lambertian scatterer has $N_{\mathrm{mix}} = 1$. The largest eigenvalue of $B_R$ describes the light after many passes. This light escapes with a probability of $p_\infty = 1-\rho(B_R)$ per encounter, where $\rho$ is the spectral radius [23,24]. We define the late-pass escape ratio as $\xi = p_\infty/\tilde p$, where $\tilde p$ is the corresponding probability for light with uniform radiance. This ratio equals unity when long-lived light escapes as easily as Lambertian light and tends to zero when some angles neither escape nor mix. The spectral quantities are computed on the bins below 85°, as described in Note S1 of Supplement 1.

**Fig. 1.** Multi-pass extraction as a Markov chain on angle. (a) Schematic of the model. (b)–(d) Return matrices $M$ of (b) a planar interface, (c) the hemispherical MLA, and (d) a scattering layer. (e) Angular distribution of the light that meets the hemispherical MLA at encounter $j$. (f) Escape probability at encounter $j$ relative to $p$, where the open symbols show $p_{\mathrm{mix}}/p$ and the dotted line shows the FR model. *(`fig1_mixing_concept.png`)*

### 2.2 Reduced model and complete mixing at the escape interface

The reduced model keeps the $p(\theta)$ of a structure and replaces $M$ with a one-parameter family. A fraction $1-m$ of the returned light keeps its angle, and a fraction $m$ is spread over a fixed distribution $q$ as

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right]. \tag{5}
$$

Equation (3) holds only for $q \propto w(1-p)$; the light spread into an angle is then proportional to the light returned by this angle. The same form appears in reciprocal reflectance models [25], whereas earlier multi-pass extraction models spread the returned light isotropically [11,15], which violates Eq. (3) whenever $p$ depends on angle. Because the kernel of Eq. (5) has the eigenvalues 1 and $1-m$, we set $m = 1-\lambda_2$.

For $m = 1$, the structure mixes the returned light completely, and this light escapes at the next encounter with a probability of

$$
p_{\mathrm{mix}} = \sum_r q_r p_r = p - \frac{\sigma_p^2}{1-p}, \qquad \sigma_p^2 = \frac{\sum_k w_k (p_k - p)^2}{\sum_k w_k}, \tag{6}
$$

while the multi-pass series takes the closed form

$$
\eta_{\mathrm{mix}} = \mathbf{p}\cdot\mathbf{v}^{(0)} + \frac{(\mathbf{p}\cdot\mathbf{u})\,[(\mathbf{1}-\mathbf{p})\cdot\mathbf{v}^{(0)}]}{1-(\mathbf{1}-\mathbf{p})\cdot\mathbf{u}}, \qquad u_k = R_k\,q_k . \tag{7}
$$

Equation (7) equals the FR model only if $p(\theta)$ is constant. The variance term is large for refractive surfaces because they transmit inside the critical cone and reflect everything outside it. For example, a planar interface has $p$ = 0.40 but $p_{\mathrm{mix}}$ = 0.06. Since the FR model assigns the probability $p$ to the mixed light, it overestimates the gain from mixing. Note that Eq. (7) is a reference value for rank-one reciprocal mixing rather than an upper bound.

### 2.3 Mixing in transmission and at the device

The penalty in Eq. (6) arises because the escape interface itself does the mixing. Reciprocity weights the redistribution in reflection by the return probability of the reflecting element, which is $1-p$ in this case. Mixing elsewhere is weighted differently, as summarized in Table 1. A lossless, non-reflecting layer that mixes the light in transmission spreads it in proportion to $w$ at every pass. The mixed light then escapes with the probability $p$, and the FR model is recovered exactly, regardless of the position of the layer in an index-matched substrate.

The device can also mix the light in reflection, for example at a diffusely reflecting electrode. In this case, the device reflection becomes a kernel $D$ with column sums of $R_i$ and with $w_i D(i\to j) = w_j D(j\to i)$. The rank-one kernel that meets both conditions is

$$
D(i\to j) = \frac{R_i\,R_j\,w_j}{\sum_k w_k R_k} . \tag{8}
$$

This kernel spreads the light in proportion to $wR$, and the mixed light escapes at the next encounter with a probability of

$$
p_D = \frac{\sum_k w_k R_k\,p_k}{\sum_k w_k R_k} = p + \frac{\mathrm{Cov}_w(R,p)}{\langle R\rangle_w}, \tag{9}
$$

where the covariance and the mean are étendue weighted. No variance term appears, and $p_D = p$ if $R$ does not depend on angle. However, if the device absorbs most strongly inside the escape cone, as strong microcavities do, the covariance is negative and mixing at the device is penalized as well (Note S4).

**Table 1.** Escape probability of completely mixed light for three placements of the mixing, and the resulting $\eta_{\mathrm{sta}}$ of a planar interface in the Al device at 550 nm

| Where the light is mixed | Light is spread in proportion to | Escape probability of the mixed light | $\eta_{\mathrm{sta}}$, planar interface |
|---|---|---|---|
| Escape interface, in reflection | $w(1-p)$ | $p - \sigma_p^2/(1-p)$ | 0.488 |
| Lossless volume, in transmission | $w$ | $p$ | 0.798 |
| Device, in reflection | $wR$ | $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$ | 0.770 |

### 2.4 Optical simulation

The dipole emission and the stack reflectance were calculated with the Chance–Prock–Silbey model extended to birefringent layers and with the transfer-matrix method, using measured optical constants and an incoherent substrate [26–28]. The reference stack was Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / glass with a refractive index of $n$ = 1.51. A second stack replaced Al with Ag, and a low-loss stack followed Ref. [13]. The $A'$ values of the three stacks at 550 nm were 0.167, 0.044, and 0.026, respectively. The results are given at 550 nm or are spectrally weighted over the visible range with the emitter spectrum.

The library consisted of 38 structures that were index matched to the substrate, namely a planar interface, six non-absorbing Henyey–Greenstein (HG) scattering layers [29], nine hexagonal MLAs with hemispherical or freeform lenslets, 17 random MLAs with independent lenslet profiles, heights, or positions, and five Gaussian wavy surfaces (Note S2). The MLAs and textures were traced with a height-field Monte Carlo (MC) ray tracer, and the escape probability of the hemispherical MLA agreed with the LightTools result. Ray optics is adequate for lenslets with a radius of 10 µm [30]. The HG layers were traced by MC and combined with the planar exit face by the adding method [31]. The bin width, the number of rays, and the residual asymmetry of the traced matrices changed $\eta_{\mathrm{sta}}$ by at most 0.4 pp, and a separate treatment of the two polarizations changed the gain from complete mixing by at most 0.4 pp.

### 2.5 Device fabrication and characterization

Orange and green bottom-emitting OLEDs with a circular pixel of 1 mm radius were fabricated on glass substrates with a size of 25.4 mm and a thickness of 0.7 mm. The orange stack was ITO (150 nm) / TAPC (50 nm) / TCTA (10 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / B3PyMPM (45 nm) / LiF (1 nm) / Al (100 nm), and the green devices used Ir(ppy)$_2$acac in the same host. For the MLA devices, a second glass of the same thickness was attached with index-matching oil. This glass carried a hexagonal array of hemispherical or freeform lenslets on a 3 × 3 mm patch centered on the pixel, and the substrate edges were blackened.

The angular intensity $I(\theta)$ was measured with a photodiode goniometer from 0° to 80° in steps of 10°. The total EQE was obtained by multiplying the on-axis EQE of a Lambertian emitter by the angular factor $F = 2\int_0^{\pi/2}[I(\theta)/I(0)]\sin\theta\,d\theta$, integrated with a spline through the measured points.

The devices were modeled by extending the chain to the lateral position (Note S6). In this model, the light moves sideways by $2t\tan\theta$ between two encounters with the front face of a substrate of thickness $t$, the back reflectance follows the device layout, and the light that reaches the edges is lost. The cavity and one internal efficiency per emitter were calibrated on the bare devices. The MLA devices were then predicted without further parameters, and their angular factors were integrated in the same way as the measured ones. The model values are given as ranges over the orientation of the lens lattice and the goniometer geometry.

---

## 3. Results and discussion

### 3.1 Reflective and transmissive mixing

Figure 1(e) shows the angular distribution of the light that meets the hemispherical MLA at successive encounters in the Al device at 550 nm. While 28% of the light lies beyond 60° at the first encounter, this share increases to 83% after ten returns, and the escape probability per encounter falls to a third of $p$, as shown in Fig. 1(f). Repeated passes thus concentrate the light where it cannot escape. Structures that mix within one or two encounters do not return to $p$ either, and an HG layer and a steep wavy surface settled at their $p_{\mathrm{mix}}$, as predicted by Eq. (6).

Figure 2(a) compares $p_{\mathrm{mix}}$ with $p$ for all structures. The ratio $p_{\mathrm{mix}}/p$ ranges from 0.81 to 0.93 for the MLAs and falls to 0.14 for the planar interface, which indicates that the penalty is largest for a surface relying on TIR. Mixing at the device avoids this penalty in the three stacks studied, because their reflectance is nearly uncorrelated with $p(\theta)$. As shown in Fig. 2(b), the advantage of mixing at the device over mixing at the escape interface therefore increases with $p - p_{\mathrm{mix}}$. For the planar interface, complete mixing raised $\eta_{\mathrm{sta}}$ from 0.384 to only 0.488 at the escape interface, compared with 0.770 at the device and 0.798 in transmission, whereas the advantage amounts to a few pp for the MLAs.

Figure 2(e) shows the gain of each structure for complete mixing of the light with its own $p(\theta)$ unchanged. For the 37 non-planar structures, Eq. (7) changed $\eta_{\mathrm{sta}}$ by −0.4 to +1.8 pp in the Al device at 550 nm and by at most 2.2 pp with spectral weighting in any of the three devices. In contrast, the FR model predicts gains of up to 19 pp for the same structures. Reflective mixing is more helpful when the cavity emits mainly into angles with poor transmission through the structure. At 500 nm, where the first-pass escape probability of the hemispherical MLA is well below $p$, complete mixing would add 6 to 9 pp, and the random MLAs and wavy surfaces realize most of this gain (Table S9).

Figure 2(d) shows that transmissive mixing does not require an ideal randomizer. A forward-scattering HG volume with an asymmetry parameter of $g$ ≥ 0.8 and a reduced optical thickness of $\tau^*$ ≈ 0.5 under the planar exit face reached $\eta_{\mathrm{sta}}$ = 0.764 in the Al device. This value is 92% of the way from the planar interface to the ideal transmissive value, and the same volume recovered 97% and 99% in the Ag and low-loss devices, respectively. However, thicker HG volumes lose this advantage because they reflect more light back to the OLED.

This distinction is not specific to OLEDs. Figure 2(c) shows the LEE of a slab with an internal isotropic source and a back reflectance of 0.90. At $n$ = 1.5, complete mixing at a planar exit face raised the LEE from 0.44 to 0.60, whereas mixing at the back reflector, equivalent to transmissive mixing in this case, gave 0.87. This ordering persists up to $n$ = 3.5. Textured exit faces have a small $\sigma_p^2$ and therefore come close to the FR result without any mixing. The same ordering also appears in light trapping, the reciprocal problem. In a weakly absorbing slab, mixing only the light that the front face reflects never gave a longer path than a Lambertian back reflector (Tables S11 and S12).

**Fig. 2.** Dependence of the LEE on where the light is mixed. (a) $p_{\mathrm{mix}}$ against $p$ for all structures. (b) Difference in $\eta_{\mathrm{sta}}$ between complete mixing at the device and at the escape interface. (c) LEE of a slab with a back reflector against its refractive index. (d) $\eta_{\mathrm{sta}}$ of HG volumes under a planar exit face. (e) Gain in $\eta_{\mathrm{sta}}$ for complete mixing by the structure (filled symbols) and in the FR model (open symbols). Panels (a), (b), (d), and (e) are for the Al device at 550 nm. *(`fig2_complete_mixing.png`)*

### 3.2 Reduced model and structure map

Figure 3(a) compares the reduced model with the full series of Eq. (2) for 342 cases, comprising 38 structures, three devices, and three wavelengths. With $m = 1-\lambda_2$ and reciprocal spreading, the mean absolute error was 0.30 pp and the largest error was 3.0 pp, without any concentration in a particular device or wavelength, as shown in Fig. 3(b). For comparison, the FR model erred by up to 60 pp for the planar interface, and a reduced model with isotropic spreading erred by up to 10.5 pp for shallow wavy surfaces. In addition, direct ray tracing in the polar and azimuthal angles agreed with the full series to within 0.6 pp (Note S3).

Figure 3(c) shows the limits of the reduced model, which we tested on 21 additional structures chosen after the model had been fixed. The model fails when the chain is reducible, i.e., when a set of angles is closed under reflection and $\lambda_2$ exceeds 0.999. This is the case for cones, pyramids, and prisms with a side angle of 45°, where near-normal light is retroreflected and the errors were about 20 pp. For the pyramids and prisms, this failure was predicted from $M$ before the series was evaluated. A second diagnostic, comparing $\lambda_2$ with the relaxation rate weighted by the first returned light, $\lambda_{\mathrm{eff}}$, identified the cones with a side angle of 60° but was not reliable in the prospective test. For the remaining faceted structures, the error reached 7.5 pp, whereas the full matrix has none of these restrictions.

Figure 4(a) shows the structures in the $(p, N_{\mathrm{mix}})$ plane, where the structure classes occupy separate regions. The MLAs have the highest escape probabilities of 0.42 to 0.43 with $N_{\mathrm{mix}}$ values of 3 to 8, whereas the HG layers relax within about two encounters but have lower escape probabilities. Randomizing the lenslet profiles of an MLA lowered $N_{\mathrm{mix}}$ from 7.0 to 1.9 without lowering $p$, as shown in Fig. 4(b), and $\xi$ accordingly increased from 0.39 to 0.80, as shown in Fig. 4(c). However, $\eta_{\mathrm{sta}}$ changed by less than 1 pp, because the gain from faster mixing was offset by the accompanying change of $p(\theta)$. The efficiency is therefore determined mainly by $p$, as indicated by the nearly vertical contours in Fig. 4(a).

Figure 4(a) also explains why optimizing the lens shape of an MLA gives only small gains. A hemispherical MLA already has an escape probability close to the étendue limit for uniform radiance, $1/n^2$ = 0.439. Although the non-uniform emission of the cavity would allow a first-pass value of up to 0.540 [32], refractive lenslets fill each outgoing direction from a wide range of internal angles and thus transmit about 0.43 in either case. Consistently, in a LightTools benchmark on the device platform of Ref. [8], optimized freeform lenslets exceeded the optimized hemisphere by only 1.3% in total EQE (Note S5). Approaching the bound would require a surface with angle-selective transmission.

**Fig. 3.** Accuracy of the reciprocal reduced model. (a) $\eta_{\mathrm{sta}}$ from the reduced model against the full series for 342 cases. (b) Error in (a) for each device (LL, low loss) and wavelength, where the bars give the mean absolute error. (c) Largest error against $\lambda_2 - \lambda_{\mathrm{eff}}$ for the library structures and the additional structures, where the crosses mark reducible chains. *(`fig3_reduced_model.png`)*

**Fig. 4.** Escape probability and mixing number of the structures. (a) $p$ and $N_{\mathrm{mix}}$ of all non-planar structures, with the contours of $\eta_{\mathrm{sta}}$ from the reduced model for the $p(\theta)$ shape of the hemispherical MLA (Al device, 550 nm). (b) $N_{\mathrm{mix}}$ against the normalized randomization parameters. (c) Late-pass escape ratio $\xi$ against $N_{\mathrm{mix}}$. *(`fig4_structure_map.png`)*

### 3.3 Fabricated devices

Figure 5(a) shows the layout of the fabricated devices. The MLA patch is only 3 mm wide, whereas the substrate with the attached lens glass is 1.4 mm thick. As a result, most of the light returned by the MLA lands outside the patch, and the MLA acts mainly on the first-pass light. The experiment therefore tests how the model redistributes this light in angle.

Figures 5(b)–5(d) and Table 2 compare the measured angular behavior with the model. The model predicted the bare, hemispherical, and freeform devices of both emitters without any adjustment. The hemispherical MLA concentrated the emission toward the normal and lowered the angular factor to about 0.7, whereas the freeform MLAs spread the first-pass light and kept the factor close to unity. All measured means lie inside the predicted ranges, except for the orange freeform devices, whose mean is 0.03 above the range. The mean angular profiles deviate from the centers of the predicted bands by at most 0.04 in rms.

**Table 2.** Angular factors of the fabricated devices, where the model values are ranges over the measurement geometry and the measured values are means ± standard deviations

| Emitter | Device | Model | Measured | Sample |
|---|---|---|---|---|
| Orange | Bare | 0.99–1.02 | 1.01 ± 0.05 | 5 devices |
| Orange | Hemispherical MLA | 0.62–0.68 | 0.66 ± 0.04 | 5 devices |
| Orange | Freeform MLA | 0.92–1.04 | 1.07 ± 0.07 | 4 devices |
| Green | Bare | 1.00–1.02 | 1.02 | 2 devices |
| Green | Hemispherical MLA | 0.65–0.71 | 0.71 ± 0.02 | 1 device, 6 scans |
| Green | Freeform MLA | 0.83–0.91 | 0.87 ± 0.04 | 3 devices |

Figure 5(e) compares the on-axis and total EQEs of the MLA devices with those of the bare devices. Because the hemispherical MLA narrows the emission, its on-axis gain overstates the LEE, and while the on-axis EQE was 1.2 times that of the bare devices, the total EQE was only 0.81 to 0.85 times. The measured EQEs of the MLA devices were 84% to 91% of the lossless prediction. A single angle-independent transmission of 0.925 per crossing of the oil-coupled lens glass, fitted to the orange hemispherical devices, reproduced this offset for the other three device types.

The model also indicates the required width of an MLA patch, because the patch has to extend several substrate thicknesses beyond the pixel in order to extract the recycled light. A 10 mm patch gives 80% of the gain of a 20 mm patch, whereas a 3 mm patch lowers the total EQE, which agrees with earlier observations on laboratory samples and pixelated displays [8,17]. Therefore, the angle-integrated EQE should be reported together with the on-axis value.

**Fig. 5.** Fabricated devices and finite-geometry model. (a) Substrate layout with the pixel and the MLA patch. (b), (c) Normalized angular intensity of the (b) orange and (c) green devices, where the bands show the model and the symbols show the measurements. (d) Angular factor. (e) Ratio of the MLA devices to the bare devices for the on-axis and total EQE, where the dashed outlines show the lossless model and the filled boxes show the model with the fitted lens-glass transmission. *(`fig5_experiment.png`)*

### 3.4 Design implications

These results suggest the following guidelines for the design of outcoupling structures. Before a structure is optimized for mixing, Eq. (6) can be used to estimate the possible gain from mixing, about 2 pp at most for the MLAs studied here. Larger gains require mixing at a position without the reciprocity penalty, i.e., in a forward-scattering volume or at a device whose reflectance is uncorrelated with $p(\theta)$. The latter adds 1.5 to 4.1 pp for the MLAs in the Al device. Furthermore, structures should be ranked by $p$ first, because $N_{\mathrm{mix}}$ and $\xi$ matter mainly when the first-pass distribution does not match $p(\theta)$. Although random MLAs and steep wavy surfaces combine a high $p$ with fast mixing, their spectrally weighted advantage over periodic MLAs is at most 2 pp. They remain attractive mainly because they reduce the angular color shift [33]. These conclusions also hold with spectral weighting, for exit media with refractive indices of up to 3.5, for finite substrates, and for partly diffuse electrodes (Table S13).

---

## 4. Conclusion

In this work, we showed that optical reciprocity ties the value of angular mixing to the position of the mixing. The light mixed in reflection by the escape interface escapes with a probability of $p - \sigma_p^2/(1-p)$, whereas the light mixed in transmission through a lossless volume escapes with $p$, as assumed in the FR model of statistical ray optics. For 37 scattering, microlens, and textured structures, complete mixing by the structure added at most 2.2 pp to the spectrally weighted LEE, whereas a forward-scattering volume recovered 92% to 99% of the FR gain over a planar interface. A one-parameter reciprocal model reproduced the full calculation with a mean absolute error of 0.30 pp for curved and random structures, and fabricated devices with two emitters confirmed the predicted angular redistribution of the first-pass light. It should be noted that the analysis applies to passive external structures index matched to the substrate and large enough for ray optics. We believe that these results will be helpful for designing outcoupling structures. A substantially higher LEE requires mixing in transmission or at the device, MLA patches much wider than the substrate thickness, or angle-selective structures.

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

6. W. H. Koo, S. M. Jeong, F. Araoka, et al., "Light extraction from organic light-emitting diodes enhanced by spontaneously formed buckles," Nat. Photonics **4**, 222–226 (2010).

7. M. Kovačič, "Design of external microtextures for efficient light outcoupling in OLEDs with different preferential orientation of emission dipoles," Opto-Electron. Rev. **30**, e141542 (2022).

8. M. Kim, J. Kim, and S. Yoo, "Near-planar light outcoupling structures with finite lateral dimensions for ultra-efficient and optical crosstalk-free OLED displays," Nat. Commun. **16**, 11606 (2025).

9. J. Kim, E. Kim, J. Park, et al., "Toward near-foldable surface light sources with ultimate efficiency: ultrathin substrates embedded with micron-scale inverted lens arrays," ACS Photonics **10**, 1775–1782 (2023).

10. E. Yablonovitch, "Statistical ray optics," J. Opt. Soc. Am. **72**, 899–907 (1982).

11. I. Schnitzer, E. Yablonovitch, C. Caneau, et al., "30% external quantum efficiency from surface textured, thin-film light-emitting diodes," Appl. Phys. Lett. **63**, 2174–2176 (1993).

12. M. A. Green, "Lambertian light trapping in textured solar cells and light-emitting diodes: analytical solutions," Prog. Photovoltaics **10**, 235–241 (2002).

13. J. Kim, S. Lee, S. Yoo, et al., "Minimizing round-trip loss for near-unity quantum efficiency in organic light-emitting diodes," arXiv, arXiv:2610.12065 (2026).

14. J. J. Shiang and A. R. Duggal, "Application of radiative transport theory to light extraction from organic light emitting diodes," J. Appl. Phys. **95**, 2880–2888 (2004).

15. K. Neyts and A. Ullan Nieto, "Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices," J. Opt. Soc. Am. A **23**, 1201 (2006).

16. J. Eisenlohr, N. Tucher, O. Höhn, et al., "Matrix formalism for light propagation and absorption in thick textured optical sheets," Opt. Express **23**, A502–A518 (2015).

17. M. Kovačič, D. Samigullina, F. Bouchard, et al., "Analysis and optimization of light outcoupling in OLEDs with external hierarchical textures," Opt. Express **29**, 23701 (2021).

18. H. W. Deckman, C. B. Roxlo, and E. Yablonovitch, "Maximum statistical increase of optical absorption in textured semiconductor films," Opt. Lett. **8**, 491–493 (1983).

19. F. J. J. Clarke and D. J. Parry, "Helmholtz reciprocity: its validity and application to reflectometry," Light. Res. Technol. **17**, 1–11 (1985).

20. E. Veach, "Robust Monte Carlo methods for light transport simulation," Ph.D. thesis, Stanford University (1997).

21. M. C. W. van Rossum and T. M. Nieuwenhuizen, "Multiple scattering of classical waves: microscopy, mesoscopy, and diffusion," Rev. Mod. Phys. **71**, 313–371 (1999).

22. D. A. Levin and Y. Peres, *Markov Chains and Mixing Times*, 2nd ed. (American Mathematical Society, 2017).

23. J. N. Darroch and E. Seneta, "On quasi-stationary distributions in absorbing discrete-time finite Markov chains," J. Appl. Probab. **2**, 88–100 (1965).

24. E. G. Altmann, J. S. E. Portela, and T. Tél, "Leaking chaotic systems," Rev. Mod. Phys. **85**, 869–918 (2013).

25. M. Ashikhmin and P. Shirley, "An anisotropic Phong BRDF model," J. Graph. Tools **5**, 25–32 (2000).

26. R. R. Chance, A. Prock, and R. Silbey, "Molecular fluorescence and energy transfer near interfaces," Adv. Chem. Phys. **37**, 1–65 (1978).

27. C.-K. Moon, S.-Y. Kim, J.-H. Lee, et al., "Luminescence from oriented emitting dipoles in a birefringent medium," Opt. Express **23**, A279 (2015).

28. S. Mladenovski, S. Hofmann, S. Reineke, et al., "Integrated optical model for organic light-emitting devices," J. Appl. Phys. **109**, 083114 (2011).

29. L. G. Henyey and J. L. Greenstein, "Diffuse radiation in the galaxy," Astrophys. J. **93**, 70–83 (1941).

30. K. Tang, R. A. Dimenna, and R. O. Buckius, "Regions of validity of the geometric optics approximation for angular scattering from very rough surfaces," Int. J. Heat Mass Transfer **40**, 49–59 (1997).

31. H. C. van de Hulst, *Multiple Light Scattering: Tables, Formulas, and Applications* (Academic Press, 1980).

32. J. Muschaweck and C. Wiesmann, "Extracting light out of LEDs," Adv. Opt. Technol. **2**, 291–300 (2013).

33. J. Lim, S. S. Oh, D. Y. Kim, et al., "Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens array," Opt. Express **14**, 6564 (2006).
