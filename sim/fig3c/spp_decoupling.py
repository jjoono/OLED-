#!/usr/bin/env python3
"""Why the n_e curves of Fig. 3(c)(ii) cross at d_ETL ~ 140 nm, and what that
means for driving the plasmon loss to zero.

Four tests, in order:
  1  is it the TCO?                  -> no; an index-matched lossless electrode
                                        gives the same ordering
  2  where does the residual sit?    -> in each curve's own plasmon peak, not in
                                        a TCO mode and not in the near field
  3  how fast does it decouple?      -> as exp(-2 kappa d) with kappa the field
                                        decay constant of the mode INSIDE the
                                        layer that fills the emitter-metal gap,
                                        kappa = (n_o/n_e) sqrt(k_SPP^2 - n_e^2 k0^2).
                                        Asymptotically this is ~lambda/12 for
                                        every n_e: lowering n_e lowers k_SPP by
                                        almost exactly as much as it lowers n_e,
                                        so the exponential rate does not improve.
                                        The thin-ETL gain is a finite-thickness
                                        effect (the mode is still hybridised with
                                        the isotropic layers, so its index is
                                        pulled up and kappa with it) and it fades
                                        as the layer thickens -- hence the crossing.
  4  must the low-n_e layer fill the gap?  -> yes
"""
import numpy as np, sys
import cps2, modes, nspp, materials as M

AG = (M.AG, M.AG, 100.0)
BOT = [(1.8, 1.8, 50.0), (M.ITO, M.ITO, 50.0)]


def stack(above, ito=M.ITO):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=above,
                      below=[(1.8, 1.8, 50.0), (ito, ito, 50.0)], n_sub=1.8)


def etl(ne, d, ito=M.ITO):
    return stack([(1.8, ne, d), AG], ito)


def hdr(t):
    print('\n' + t + '\n' + '-' * len(t))


hdr('1. is the crossing caused by the transparent electrode?')
print('   SPP fraction (%) at d_ETL = 150 / 200 nm')
print('   %-34s %-13s %-13s %s' % ('bottom electrode', 'n_e=1.80', '1.70', '1.60'))
for lab, ito in (('ITO 1.8636 + 0.0032i (baseline)', M.ITO),
                 ('lossless ITO 1.8636', 1.8636 + 0j),
                 ('index-matched 1.80, lossless', 1.80 + 0j)):
    c = ['%5.2f /%5.2f' % tuple(100 * cps2.solve(etl(ne, d, ito))['spp']
                                for d in (150.0, 200.0)) for ne in (1.80, 1.70, 1.60)]
    print('   %-34s %s' % (lab, ' '.join('%-13s' % x for x in c)))
print('   -> the ordering survives an index-matched lossless electrode, so no.')

hdr('2. where in k_x/k0 does the remaining loss sit?')
BINS = [(1.801, 1.85), (1.85, 1.95), (1.95, 2.15), (2.15, 5.40)]
print('   fraction of the total dissipated power, d_ETL = 200 nm')
print('   n_e  ' + '  '.join('%10s' % ('%.2f-%.2f' % b) for b in BINS) + '     total')
for ne in (1.80, 1.70, 1.60):
    S = etl(ne, 200.0)
    p = []
    for a, b in BINS:
        g = np.linspace(a, b, 4001)
        with np.errstate(all='ignore'):
            tm, te = cps2.spectrum(S, g)
        p.append(100 * np.trapezoid(tm + te, g))
    print('   %.2f  ' % ne + '  '.join('%10.3f' % x for x in p)
          + '   %7.3f' % (100 * cps2.solve(S)['spp']))
print('   -> all of it is the plasmon of that particular n_e; the band just above')
print('      the light line (where a TCO mode would sit) and the near field are flat.')

hdr('3. the decoupling rate')
print('   complex TM pole of the stack, and 1/(2 Re kappa) inside the ETL')
print('   n_e  d_ETL     n_mode              1/(2 kappa) in the ETL')
for ne, lo, hi in ((1.80, 1.95, 2.15), (1.70, 1.86, 2.00), (1.60, 1.801, 1.90)):
    for d in (60.0, 100.0, 150.0, 200.0, 250.0):
        p = modes.tm_pole(etl(ne, d), lo, hi)
        _, L = modes.kappa(p, 1.8, ne)
        print('   %.2f %5.0f    %.5f %+0.5fi        %6.1f nm' % (ne, d, p.real, p.imag, L))
    print()
print('   asymptotic value from the analytic plasmon index (semi-infinite ETL):')
print('     n_e    n_SPP     1/(2 kappa)')
for ne in (1.80, 1.75, 1.70, 1.65, 1.60):
    n = nspp.n_spp(M.AG, 1.8, ne)
    q = (1.8 / ne) ** 2 * (n ** 2 - ne ** 2)
    print('     %.2f   %.4f    %6.1f nm' % (ne, n, 550.0 / (4 * np.pi * np.sqrt(q))))
print('   -> ~lambda/12 whatever n_e is. Lowering n_e does NOT speed up the')
print('      exponential; at 60 nm it does (39 nm vs 45 nm) only because the mode')
print('      is still hybridised, and that advantage decays away with thickness.')

hdr('4. must the low-n_e material fill the emitter-metal gap?')
print('   SPP fraction (%) against the total gap between EML and Ag')
print('   gap   A: n_e=1.6 fills it   B: n_e=1.6 only 20 nm at Ag   C: isotropic 1.8')
for d in (60.0, 100.0, 150.0, 200.0, 250.0):
    a = 100 * cps2.solve(etl(1.60, d))['spp']
    b = 100 * cps2.solve(stack([(1.8, 1.8, d - 20.0), (1.8, 1.60, 20.0), AG]))['spp']
    c = 100 * cps2.solve(etl(1.80, d))['spp']
    print('   %4.0f          %6.2f                  %6.2f                    %6.2f'
          % (d, a, b, c))
print('   -> yes. A low-n_e skin at the metal is worse than a uniform low-n_e layer')
print('      everywhere, and worse than plain isotropic once the gap exceeds ~150 nm.')

hdr('5. the only route to zero')
print('   SPP fraction (%) across the threshold, n_o = 1.80')
print('   n_e   n_SPP     d=100   d=150   d=200   d=250')
A = np.genfromtxt('fig3c_sweep.csv', delimiter=',', names=True, dtype=None,
                  encoding='utf-8')
F = A[A['series'] == 'no1.80']
for ne in (1.70, 1.64, 1.60, 1.58, 1.56, 1.50):
    m = F[np.isclose(F['n_e'], ne)]
    o = np.argsort(m['d_ETL_nm'])
    v = [100 * np.interp(d, m['d_ETL_nm'][o], m['spp'][o])
         for d in (100, 150, 200, 250)]
    print('   %.2f  %.4f  ' % (ne, nspp.n_spp(M.AG, 1.8, ne))
          + ' '.join('%7.2f' % x for x in v))
print('   -> thickness alone buys exp(-d / 46 nm) whatever n_e is. The loss only')
print('      collapses when n_SPP drops below the organic/substrate light line and')
print('      the mode stops being bound at all.')
