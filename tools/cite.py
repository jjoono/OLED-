"""Reference database (Optica style) and citation resolver for the manuscript sources.

Usage:  python3 tools/cite.py <source.md> <output.md>

The source uses keys, e.g. [@greiner2007;@kim2026], which are numbered in order of first appearance.
Consecutive numbers are compressed as [1,2] (two) or [1–3] (three or more).  A block '## References'
with the cited entries is appended.  Every entry was checked against Crossref / OpenAlex / arXiv
(see review/reference_check.md); 'doi' is kept here for that record and is not printed in the list.
"""
import re, sys

REFS = {
 'brutting2013': ('W. Brütting, J. Frischeisen, T. D. Schmidt, et al., "Device efficiency of organic light-emitting diodes: progress by improved light outcoupling," Phys. Status Solidi A **210**, 44–65 (2013).', '10.1002/pssa.201228320'),
 'greiner2007': ('H. Greiner, "Light extraction from organic light emitting diode substrates: simulation and experiment," Jpn. J. Appl. Phys. **46**, 4125–4137 (2007).', '10.1143/JJAP.46.4125'),
 'moller2002': ('S. Möller and S. R. Forrest, "Improved light out-coupling in organic light emitting diodes employing ordered microlens arrays," J. Appl. Phys. **91**, 3324–3327 (2002).', '10.1063/1.1435422'),
 'qu2018': ('Y. Qu, J. Kim, C. Coburn, et al., "Efficient, nonintrusive outcoupling in organic light emitting devices using embedded microlens arrays," ACS Photonics **5**, 2453–2458 (2018).', '10.1021/acsphotonics.8b00255'),
 'song2018': ('J. Song, K.-H. Kim, E. Kim, et al., "Lensfree OLEDs with over 50% external quantum efficiency via external scattering and horizontally oriented emitters," Nat. Commun. **9**, 3207 (2018).', '10.1038/s41467-018-05671-x'),
 'levell2013': ('J. W. Levell, S. Harkema, R. K. Pendyala, et al., "Device reflectivity as a simple rule for predicting the suitability of scattering foils for improved OLED light extraction," Proc. SPIE **8829**, 88291L (2013).', '10.1117/12.2023021'),
 'kovacic2022': ('M. Kovačič, "Design of external microtextures for efficient light outcoupling in OLEDs with different preferential orientation of emission dipoles," Opto-Electron. Rev. **30**, e141542 (2022).', '10.24425/opelre.2022.141542'),
 'sun2006': ('Y. Sun and S. R. Forrest, "Organic light emitting devices with enhanced outcoupling via microlenses fabricated by imprint lithography," J. Appl. Phys. **100**, 073106 (2006).', '10.1063/1.2356904'),
 'kim2021': ('S. Kim, J. M. Shin, J. Lee, et al., "Inverse design of organic light-emitting diode structure based on deep neural networks," Nanophotonics **10**, 4533–4541 (2021).', '10.1515/nanoph-2021-0434'),
 'zhao2026': ('Z.-R. Zhao, K.-P. Song, H. Liu, et al., "Deep Q-network-guided TiO$_2$/NOA61 microstructures for light outcoupling and angular emission redistribution in blue top-emitting OLEDs," Adv. Funct. Mater. **36**, e78939 (2026).', '10.1002/adfm.78939'),
 'kim2025': ('M. Kim, J. Kim, and S. Yoo, "Near-planar light outcoupling structures with finite lateral dimensions for ultra-efficient and optical crosstalk-free OLED displays," Nat. Commun. **16**, 11606 (2025).', '10.1038/s41467-025-66538-6'),
 'kim2023': ('J. Kim, E. Kim, J. Park, et al., "Toward near-foldable surface light sources with ultimate efficiency: ultrathin substrates embedded with micron-scale inverted lens arrays," ACS Photonics **10**, 1775–1782 (2023).', '10.1021/acsphotonics.3c00017'),
 'li2019': ('W. Li, Y.-Q. Li, Y. Shen, et al., "Releasing the trapped light for efficient silver nanowires-based white flexible organic light-emitting diodes," Adv. Opt. Mater. **7**, 1900985 (2019).', '10.1002/adom.201900985'),
 'kim2026': ('J. Kim, S. Lee, S. Yoo, et al., "Minimizing round-trip loss for near-unity quantum efficiency in organic light-emitting diodes," arXiv, arXiv:2610.12065 (2026).', 'arXiv:2610.12065'),
 'yablonovitch1982': ('E. Yablonovitch, "Statistical ray optics," J. Opt. Soc. Am. **72**, 899–907 (1982).', '10.1364/JOSA.72.000899'),
 'schnitzer1993': ('I. Schnitzer, E. Yablonovitch, C. Caneau, et al., "30% external quantum efficiency from surface textured, thin-film light-emitting diodes," Appl. Phys. Lett. **63**, 2174–2176 (1993).', '10.1063/1.110575'),
 'green2002': ('M. A. Green, "Lambertian light trapping in textured solar cells and light-emitting diodes: analytical solutions," Prog. Photovoltaics **10**, 235–241 (2002).', '10.1002/pip.404'),
 'shiang2004': ('J. J. Shiang and A. R. Duggal, "Application of radiative transport theory to light extraction from organic light emitting diodes," J. Appl. Phys. **95**, 2880–2888 (2004).', '10.1063/1.1644037'),
 'eisenlohr2015': ('J. Eisenlohr, N. Tucher, O. Höhn, et al., "Matrix formalism for light propagation and absorption in thick textured optical sheets," Opt. Express **23**, A502–A518 (2015).', '10.1364/OE.23.00A502'),
 'david2013': ('A. David, "Surface-roughened light-emitting diodes: an accurate model," J. Disp. Technol. **9**, 301–316 (2013).', '10.1109/JDT.2013.2240373'),
 'altazin2015': ('S. Altazin, C. Reynaud, U. M. Mayer, et al., "Simulations, measurements, and optimization of OLEDs with scattering layer," SID Symp. Dig. Tech. Pap. **46**, 564–567 (2015).', '10.1002/sdtp.10478'),
 'kovacic2021': ('M. Kovačič, D. Samigullina, F. Bouchard, et al., "Analysis and optimization of light outcoupling in OLEDs with external hierarchical textures," Opt. Express **29**, 23701 (2021).', '10.1364/OE.428021'),
 'deckman1983': ('H. W. Deckman, C. B. Roxlo, and E. Yablonovitch, "Maximum statistical increase of optical absorption in textured semiconductor films," Opt. Lett. **8**, 491–493 (1983).', '10.1364/OL.8.000491'),
 'manzoor2020': ('S. Manzoor, M. Filipič, A. Onno, et al., "Visualizing light trapping within textured silicon solar cells," J. Appl. Phys. **127**, 063104 (2020).', '10.1063/1.5131173'),
 'battaglia2012': ('C. Battaglia, M. Boccard, F.-J. Haug, et al., "Light trapping in solar cells: when does a Lambertian scatterer scatter Lambertianly?" J. Appl. Phys. **112**, 094504 (2012).', '10.1063/1.4761988'),
 'vanrossum1999': ('M. C. W. van Rossum and T. M. Nieuwenhuizen, "Multiple scattering of classical waves: microscopy, mesoscopy, and diffusion," Rev. Mod. Phys. **71**, 313–371 (1999).', '10.1103/RevModPhys.71.313'),
 'levin2017': ('D. A. Levin and Y. Peres, *Markov Chains and Mixing Times*, 2nd ed. (American Mathematical Society, 2017).', '10.1090/mbk/107'),
 'darroch1965': ('J. N. Darroch and E. Seneta, "On quasi-stationary distributions in absorbing discrete-time finite Markov chains," J. Appl. Probab. **2**, 88–100 (1965).', '10.2307/3211876'),
 'altmann2013': ('E. G. Altmann, J. S. E. Portela, and T. Tél, "Leaking chaotic systems," Rev. Mod. Phys. **85**, 869–918 (2013).', '10.1103/RevModPhys.85.869'),
 'ryu2006': ('J.-W. Ryu, S.-Y. Lee, C.-M. Kim, et al., "Survival probability time distribution in dielectric cavities," Phys. Rev. E **73**, 036207 (2006).', '10.1103/PhysRevE.73.036207'),
 'clauss2019': ('K. Clauß, E. G. Altmann, A. Bäcker, et al., "Structure of resonance eigenfunctions for chaotic systems with partial escape," Phys. Rev. E **100**, 052205 (2019).', '10.1103/PhysRevE.100.052205'),
 'veach1997': ('E. Veach, "Robust Monte Carlo methods for light transport simulation," Ph.D. thesis, Stanford University (1997).', 'thesis (no DOI)'),
 'soler2026': ('C. Soler and K. Subr, "Spectral theory of light transport operators," ACM Trans. Graph. **45**, 1–18 (2026).', '10.1145/3774756'),
 'koenderink1983': ('J. J. Koenderink and A. J. van Doorn, "Geometrical modes as a general method to treat diffuse interreflections in radiometry," J. Opt. Soc. Am. **73**, 843 (1983).', '10.1364/JOSA.73.000843'),
 'ashikhmin2000': ('M. Ashikhmin and P. Shirley, "An anisotropic Phong BRDF model," J. Graph. Tools **5**, 25–32 (2000).', '10.1080/10867651.2000.10487522'),
 'kelemen2001': ('C. Kelemen and L. Szirmay-Kalos, "A microfacet based coupled specular-matte BRDF model with importance sampling," in *Eurographics 2001 – Short Presentations* (Eurographics Association, 2001).', '10.2312/egs.20011003'),
 'muschaweck2013': ('J. Muschaweck and C. Wiesmann, "Extracting light out of LEDs," Adv. Opt. Technol. **2**, 291–300 (2013).', '10.1515/aot-2013-0027'),
 'lenef2018': ('A. Lenef, A. Piquette, and J. Kelso, "Thermodynamics of light extraction from luminescent materials," ECS J. Solid State Sci. Technol. **7**, R3211–R3226 (2018).', '10.1149/2.0271801jss'),
 'clarke1985': ('F. J. J. Clarke and D. J. Parry, "Helmholtz reciprocity: its validity and application to reflectometry," Light. Res. Technol. **17**, 1–11 (1985).', '10.1177/14771535850170010301'),
 'neyts2006': ('K. Neyts and A. Ullan Nieto, "Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices," J. Opt. Soc. Am. A **23**, 1201 (2006).', '10.1364/JOSAA.23.001201'),
 'horn2012': ('R. A. Horn and C. R. Johnson, *Matrix Analysis*, 2nd ed. (Cambridge University Press, 2012).', '10.1017/CBO9781139020411'),
 'chance1978': ('R. R. Chance, A. Prock, and R. Silbey, "Molecular fluorescence and energy transfer near interfaces," Adv. Chem. Phys. **37**, 1–65 (1978).', '10.1002/9780470142561.ch1'),
 'wasey2000': ('J. A. E. Wasey, A. Safonov, I. D. W. Samuel, et al., "Effects of dipole orientation and birefringence on the optical emission from thin films," Opt. Commun. **183**, 109–121 (2000).', '10.1016/S0030-4018(00)00806-3'),
 'moon2015': ('C.-K. Moon, S.-Y. Kim, J.-H. Lee, et al., "Luminescence from oriented emitting dipoles in a birefringent medium," Opt. Express **23**, A279 (2015).', '10.1364/OE.23.00A279'),
 'mladenovski2011': ('S. Mladenovski, S. Hofmann, S. Reineke, et al., "Integrated optical model for organic light-emitting devices," J. Appl. Phys. **109**, 083114 (2011).', '10.1063/1.3576114'),
 'henyey1941': ('L. G. Henyey and J. L. Greenstein, "Diffuse radiation in the galaxy," Astrophys. J. **93**, 70–83 (1941).', '10.1086/144246'),
 'vandehulst1980': ('H. C. van de Hulst, *Multiple Light Scattering: Tables, Formulas, and Applications* (Academic Press, 1980).', 'book (chapter DOI 10.1016/B978-0-12-710701-1.50013-6)'),
 'koo2010': ('W. H. Koo, S. M. Jeong, F. Araoka, et al., "Light extraction from organic light-emitting diodes enhanced by spontaneously formed buckles," Nat. Photonics **4**, 222–226 (2010).', '10.1038/nphoton.2010.7'),
 'tang1997': ('K. Tang, R. A. Dimenna, and R. O. Buckius, "Regions of validity of the geometric optics approximation for angular scattering from very rough surfaces," Int. J. Heat Mass Transfer **40**, 49–59 (1997).', '10.1016/S0017-9310(96)00073-7'),
 'lim2006': ('J. Lim, S. S. Oh, D. Y. Kim, et al., "Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens array," Opt. Express **14**, 6564 (2006).', '10.1364/OE.14.006564'),
}


def _fmt(nums):
    nums = sorted(set(nums)); out = []; i = 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        if j - i >= 2: out.append(f'{nums[i]}–{nums[j]}')
        else: out += [str(n) for n in nums[i:j + 1]]
        i = j + 1
    return '[' + ','.join(out) + ']'


def resolve(text, order=None):
    """Replace [@key;@key] by numbers in order of first appearance; return (text, ordered keys)."""
    order = [] if order is None else order

    def sub(m):
        keys = [k.strip().lstrip('@') for k in m.group(1).split(';')]
        for k in keys:
            if k not in REFS: raise KeyError(k)
            if k not in order: order.append(k)
        return _fmt([order.index(k) + 1 for k in keys])
    return re.sub(r'\[((?:@[\w]+)(?:;\s*@[\w]+)*)\]', sub, text), order


def reference_list(order):
    return '\n\n'.join(f'{i + 1}. {REFS[k][0]}' for i, k in enumerate(order))


if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    body, order = resolve(open(src, encoding='utf-8').read())
    open(dst, 'w', encoding='utf-8').write(body.rstrip() + '\n\n---\n\n## References\n\n' + reference_list(order) + '\n')
    unused = [k for k in REFS if k not in order]
    print(f'{len(order)} references cited; not cited: {unused}')
