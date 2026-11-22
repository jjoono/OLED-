# 참고문헌 검증 기록 (Optics Express 투고본)

검증일: 2026-11-22. 본문 참고문헌 50건 전부를 Crossref·OpenAlex·arXiv 기록과 대조했다(저자, 제목, 학술지, 권, 쪽, 연도). 존재하지 않는 문헌은 없었다.
목록 형식은 Optica 양식(저자 3인까지 표기 후 et al., "제목," 학술지 **권**, 쪽 (연도))이고, 번호는 본문 첫 인용 순서다. 목록과 번호는 `tools/cite.py`가 `tools/manuscript_en.src.md`에서 생성한다.

## 이번 개정에서 바꾼 것

- 삭제 2건: Bortz 외 (Proc. SPIE 3139, 1997)는 Crossref와 ADS가 편집자(R. Winston)를 저자로 함께 기록하고 있어 저자 목록을 확정할 수 없었다. Ries 외 (J. Opt. Soc. Am. A 14, 1997)는 회전 대칭 광학계의 skewness 한계로, 주기 배열에 대한 본문 주장과 직접 맞지 않았다. 두 문헌이 받치던 étendue 한계 서술은 Muschaweck–Wiesmann, Lenef 외로 유지된다.
- 교체 2건: Langer (1999) → Horn & Johnson, *Matrix Analysis* (2012). Thorsos (1988) → Tang, Dimenna, Buckius (1997).
- 서술 수정: 아래 표의 비고 참조.
- 확인이 더 필요한 것: [14]는 arXiv 사전공개본이므로 투고 직전에 출판 여부를 다시 확인한다.

## 목록

| 번호 | 문헌 | DOI / 식별자 | 비고 |
|---|---|---|---|
| 1 | W. Brütting, J. Frischeisen, T. D. Schmidt, et al., "Device efficiency of organic light-emitting diodes: progress by improved light outco… | 10.1002/pssa.201228320 |  |
| 2 | H. Greiner, "Light extraction from organic light emitting diode substrates: simulation and experiment," Jpn. J. Appl. Phys. 46, 4125–4137… | 10.1143/JJAP.46.4125 |  |
| 3 | S. Möller and S. R. Forrest, "Improved light out-coupling in organic light emitting diodes employing ordered microlens arrays," J. Appl. … | 10.1063/1.1435422 |  |
| 4 | Y. Qu, J. Kim, C. Coburn, et al., "Efficient, nonintrusive outcoupling in organic light emitting devices using embedded microlens arrays,… | 10.1021/acsphotonics.8b00255 | 전극 하부 MLA(70 ± 4%). MLA 예시로만 인용하고 EQE 범위 문장에서는 제외 |
| 5 | J. Song, K.-H. Kim, E. Kim, et al., "Lensfree OLEDs with over 50% external quantum efficiency via external scattering and horizontally or… | 10.1038/s41467-018-05671-x |  |
| 6 | J. W. Levell, S. Harkema, R. K. Pendyala, et al., "Device reflectivity as a simple rule for predicting the suitability of scattering foil… | 10.1117/12.2023021 |  |
| 7 | M. Kovačič, "Design of external microtextures for efficient light outcoupling in OLEDs with different preferential orientation of emissio… | 10.24425/opelre.2022.141542 |  |
| 8 | Y. Sun and S. R. Forrest, "Organic light emitting devices with enhanced outcoupling via microlenses fabricated by imprint lithography," J… | 10.1063/1.2356904 |  |
| 9 | S. Kim, J. M. Shin, J. Lee, et al., "Inverse design of organic light-emitting diode structure based on deep neural networks," Nanophotoni… | 10.1515/nanoph-2021-0434 | 박막 스택의 DNN 역설계. "렌즈 형상 최적화"가 아니라 "스택 설계"로 서술 수정 |
| 10 | Z.-R. Zhao, K.-P. Song, H. Liu, et al., "Deep Q-network-guided TiO$_2$/NOA61 microstructures for light outcoupling and angular emission r… | 10.1002/adfm.78939 | DQN 기반 미세구조 설계. 권 36, e78939 확인 |
| 11 | M. Kim, J. Kim, and S. Yoo, "Near-planar light outcoupling structures with finite lateral dimensions for ultra-efficient and optical cros… | 10.1038/s41467-025-66538-6 | 48.0% EQE 보고. 50–60% 범위 서술과 측면 크기 논의(4.5절)에 인용 |
| 12 | J. Kim, E. Kim, J. Park, et al., "Toward near-foldable surface light sources with ultimate efficiency: ultrathin substrates embedded with… | 10.1021/acsphotonics.3c00017 |  |
| 13 | W. Li, Y.-Q. Li, Y. Shen, et al., "Releasing the trapped light for efficient silver nanowires-based white flexible organic light-emitting… | 10.1002/adom.201900985 |  |
| 14 | J. Kim, S. Lee, S. Yoo, et al., "Minimizing round-trip loss for near-unity quantum efficiency in organic light-emitting diodes," arXiv, a… | arXiv:2610.12065 | arXiv:2610.12065 (2026-10-08 등록, 저자 4인 확인). 정식 출판되면 교체 필요 |
| 15 | E. Yablonovitch, "Statistical ray optics," J. Opt. Soc. Am. 72, 899–907 (1982). | 10.1364/JOSA.72.000899 |  |
| 16 | I. Schnitzer, E. Yablonovitch, C. Caneau, et al., "30% external quantum efficiency from surface textured, thin-film light-emitting diodes… | 10.1063/1.110575 |  |
| 17 | M. A. Green, "Lambertian light trapping in textured solar cells and light-emitting diodes: analytical solutions," Prog. Photovoltaics 10,… | 10.1002/pip.404 |  |
| 18 | J. J. Shiang and A. R. Duggal, "Application of radiative transport theory to light extraction from organic light emitting diodes," J. App… | 10.1063/1.1644037 | 중간 산란 세기에서 최대라는 서술은 자체 결과(Fig. 2(d)) 근거로 바꾸고 "see also"로 인용 |
| 19 | J. Eisenlohr, N. Tucher, O. Höhn, et al., "Matrix formalism for light propagation and absorption in thick textured optical sheets," Opt. … | 10.1364/OE.23.00A502 | "상반성을 지적했다"는 표현 삭제. 재분배 행렬 사용 사례로만 인용 |
| 20 | A. David, "Surface-roughened light-emitting diodes: an accurate model," J. Disp. Technol. 9, 301–316 (2013). | 10.1109/JDT.2013.2240373 |  |
| 21 | S. Altazin, C. Reynaud, U. M. Mayer, et al., "Simulations, measurements, and optimization of OLEDs with scattering layer," SID Symp. Dig.… | 10.1002/sdtp.10478 |  |
| 22 | M. Kovačič, D. Samigullina, F. Bouchard, et al., "Analysis and optimization of light outcoupling in OLEDs with external hierarchical text… | 10.1364/OE.428021 | 실험실 크기 시료의 측면 크기·수광 면적 한계. 4.5절에 인용. 끝 쪽 번호는 미확인이라 첫 쪽만 표기 |
| 23 | H. W. Deckman, C. B. Roxlo, and E. Yablonovitch, "Maximum statistical increase of optical absorption in textured semiconductor films," Op… | 10.1364/OL.8.000491 | 흡수 증대 측정으로 서술 수정 |
| 24 | S. Manzoor, M. Filipič, A. Onno, et al., "Visualizing light trapping within textured silicon solar cells," J. Appl. Phys. 127, 063104 (20… | 10.1063/1.5131173 | 패스별 광선추적으로 서술 수정 |
| 25 | C. Battaglia, M. Boccard, F.-J. Haug, et al., "Light trapping in solar cells: when does a Lambertian scatterer scatter Lambertianly?" J. … | 10.1063/1.4761988 | Lambertian 산란 조건으로 서술 수정 |
| 26 | M. C. W. van Rossum and T. M. Nieuwenhuizen, "Multiple scattering of classical waves: microscopy, mesoscopy, and diffusion," Rev. Mod. Ph… | 10.1103/RevModPhys.71.313 |  |
| 27 | D. A. Levin and Y. Peres, Markov Chains and Mixing Times, 2nd ed. (American Mathematical Society, 2017). | 10.1090/mbk/107 |  |
| 28 | J. N. Darroch and E. Seneta, "On quasi-stationary distributions in absorbing discrete-time finite Markov chains," J. Appl. Probab. 2, 88–… | 10.2307/3211876 |  |
| 29 | E. G. Altmann, J. S. E. Portela, and T. Tél, "Leaking chaotic systems," Rev. Mod. Phys. 85, 869–918 (2013). | 10.1103/RevModPhys.85.869 |  |
| 30 | J.-W. Ryu, S.-Y. Lee, C.-M. Kim, et al., "Survival probability time distribution in dielectric cavities," Phys. Rev. E 73, 036207 (2006). | 10.1103/PhysRevE.73.036207 |  |
| 31 | K. Clauß, E. G. Altmann, A. Bäcker, et al., "Structure of resonance eigenfunctions for chaotic systems with partial escape," Phys. Rev. E… | 10.1103/PhysRevE.100.052205 |  |
| 32 | E. Veach, "Robust Monte Carlo methods for light transport simulation," Ph.D. thesis, Stanford University (1997). | thesis (no DOI) | Stanford 박사학위 논문 (DOI 없음) |
| 33 | C. Soler and K. Subr, "Spectral theory of light transport operators," ACM Trans. Graph. 45, 1–18 (2026). | 10.1145/3774756 | 실수 스펙트럼 근거로 Levin–Peres와 함께 인용. 범위(5.4절) 문단의 인용은 삭제. 2025 온라인, 2026 권호 |
| 34 | J. J. Koenderink and A. J. van Doorn, "Geometrical modes as a general method to treat diffuse interreflections in radiometry," J. Opt. So… | 10.1364/JOSA.73.000843 | 끝 쪽 번호 미확인, 첫 쪽만 표기 |
| 35 | F. J. J. Clarke and D. J. Parry, "Helmholtz reciprocity: its validity and application to reflectometry," Light. Res. Technol. 17, 1–11 (1… | 10.1177/14771535850170010301 |  |
| 36 | M. Ashikhmin and P. Shirley, "An anisotropic Phong BRDF model," J. Graph. Tools 5, 25–32 (2000). | 10.1080/10867651.2000.10487522 |  |
| 37 | C. Kelemen and L. Szirmay-Kalos, "A microfacet based coupled specular-matte BRDF model with importance sampling," in Eurographics 2001 – … | 10.2312/egs.20011003 | DOI 10.2312/egs.20011003 확인 (Eurographics 2001 Short Presentations) |
| 38 | K. Neyts and A. Ullan Nieto, "Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices," … | 10.1364/JOSAA.23.001201 | 끝 쪽 번호 미확인, 첫 쪽만 표기 |
| 39 | R. A. Horn and C. R. Johnson, Matrix Analysis, 2nd ed. (Cambridge University Press, 2012). | 10.1017/CBO9781139020411 | 신규. Rayleigh 몫·열 합 한계의 근거. 기존 Langer(1999)를 교체 |
| 40 | R. R. Chance, A. Prock, and R. Silbey, "Molecular fluorescence and energy transfer near interfaces," Adv. Chem. Phys. 37, 1–65 (1978). | 10.1002/9780470142561.ch1 |  |
| 41 | J. A. E. Wasey, A. Safonov, I. D. W. Samuel, et al., "Effects of dipole orientation and birefringence on the optical emission from thin f… | 10.1016/S0030-4018(00)00806-3 |  |
| 42 | C.-K. Moon, S.-Y. Kim, J.-H. Lee, et al., "Luminescence from oriented emitting dipoles in a birefringent medium," Opt. Express 23, A279 (… | 10.1364/OE.23.00A279 | 끝 쪽 번호 미확인, 첫 쪽만 표기 |
| 43 | S. Mladenovski, S. Hofmann, S. Reineke, et al., "Integrated optical model for organic light-emitting devices," J. Appl. Phys. 109, 083114… | 10.1063/1.3576114 |  |
| 44 | L. G. Henyey and J. L. Greenstein, "Diffuse radiation in the galaxy," Astrophys. J. 93, 70–83 (1941). | 10.1086/144246 | 저자 이니셜 L. G. (Crossref 기록의 L. C.는 오기) |
| 45 | H. C. van de Hulst, Multiple Light Scattering: Tables, Formulas, and Applications (Academic Press, 1980). | book (chapter DOI 10.1016/B978-0-12-710701-1.50013-6) | 단행본 (Academic Press) |
| 46 | W. H. Koo, S. M. Jeong, F. Araoka, et al., "Light extraction from organic light-emitting diodes enhanced by spontaneously formed buckles,… | 10.1038/nphoton.2010.7 |  |
| 47 | K. Tang, R. A. Dimenna, and R. O. Buckius, "Regions of validity of the geometric optics approximation for angular scattering from very ro… | 10.1016/S0017-9310(96)00073-7 | 신규. 기하광학 근사의 유효 범위. 기존 Thorsos(1988, 음향 Kirchhoff 근사)를 교체하고 "5 µm 이상" 문구 삭제 |
| 48 | J. Muschaweck and C. Wiesmann, "Extracting light out of LEDs," Adv. Opt. Technol. 2, 291–300 (2013). | 10.1515/aot-2013-0027 | étendue 한계의 근거 |
| 49 | A. Lenef, A. Piquette, and J. Kelso, "Thermodynamics of light extraction from luminescent materials," ECS J. Solid State Sci. Technol. 7,… | 10.1149/2.0271801jss | "이진 투과면이 한계에 도달한다"는 단정 삭제. étendue 한계의 근거로만 인용 |
| 50 | J. Lim, S. S. Oh, D. Y. Kim, et al., "Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens… | 10.1364/OE.14.006564 | 불규칙 MLA가 마이크로캐비티 소자의 각도별 색 변화를 억제한다는 내용으로 서술 수정 |
