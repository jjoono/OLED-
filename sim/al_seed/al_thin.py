"""Ultrathin Al between the bottom Ag reflector and Liq: how its absorption depends on what the Al actually becomes.
550 nm. Reflector absorption = Lambertian-weighted, unpolarised, seen from n=1.8 through Liq 1 nm.
Planar TEOLED EQE: Ag | [interlayer] | Liq 1 | Bebq2 200 | EML 25 | TAPC 80 | HATCN 10 | Ag15 | lens 1.7 MLA."""
import sys, numpy as np
sys.path[:0] = ['/home/user/OLED-/sim/isl_teoled', '/home/user/OLED-/sim/fig3c']
import isl_model as I
hw = 2.2543                                   # eV at 550 nm
AL = 0.958336812040867 + 6.68678039868499j; LIQ = 1.70; ALOX = 1.65
eAl, eAg = AL**2, I.AG**2
def drude(wp, g): return -wp**2 / (hw**2 + 1j*g*hw)
def n_of(e): n = np.sqrt(complex(e)); return n if n.imag >= 0 else -n
def al_size(d):                       # classical size effect: gamma_Al 0.13 eV -> 0.13 + hbar vF / d (1.34 eV nm)
    return n_of(eAl - drude(15.0, 0.13) + drude(15.0, 0.13 + 1.34/d))
def agal_alloy(d_al, d_mix=2.0):      # Al dissolved into the top d_mix of Ag; Linde 1.95 uOhm cm per at%
    x = 100*d_al/d_mix*0.6               # at% (Al 6.0e22 vs Ag 5.85e22 at/cm3, ~1:1; 0.6 = Al atoms per Ag site in mix, crude)
    g = 0.0174*(1.6 + 1.95*x)/1.6
    return n_of(eAg - drude(9.0, 0.0174) + drude(9.0, g)), x
TH = I.TH; w = np.cos(np.radians(TH))*np.sin(np.radians(TH)); w /= w.sum()
def A_refl(layers):                   # layers from the organic side down to the Ag substrate
    R, _ = I.rt_unpol([(LIQ, 1.0)] + layers, 1.8, I.AG, TH); return float(w @ (1 - R))
bs = np.load('/home/user/OLED-/sim/isl_teoled/mla_hemi_bsdf.npz'); BT, BR = bs['BT_1.7'], bs['BR_1.7']
tl = [(I.TAPC[0], I.TAPC[1], 80.), I.iso(I.HAT, 10.), I.iso(I.AG15, 15.)]
def eqe(inter):                        # inter: [(n, d)] between Liq and Ag, listed from Liq downwards
    bel = [I.iso(LIQ, 1.0)] + [I.iso(n, d) for n, d in inter] + [I.iso(I.AG, 100.)]
    bel = [I.iso(I.BEB, 200.)] + bel
    _, sU, _, _ = I.emission(bel, tl, 1.0, 1.7)
    Ru, _ = I.rt_unpol([(l[0], l[2]) for l in tl[::-1]] + [(I.EML[0], 25.)] + [(l[0], l[2]) for l in bel], 1.7, 1.0, TH)
    return I.solve_planar(sU, Ru, BT, BR), sU.sum()
print('reference Ag/Liq: A_refl %.2f%%, EQE %.3f (emitted into lens %.3f)' % ((100*A_refl([]),) + eqe([])))
print('%-34s %6s %7s %7s' % ('scenario (nominal Al d)', 'd', 'A_refl', 'EQE'))
for d in (0.1, 0.2, 0.3, 0.5, 1.0):
    cases = [('bulk-nk continuous Al', [(AL, d)]),
             ('size-effect Al (gamma+vF/d)', [(al_size(d), d)]),
             ('half reacted/oxidised', [(ALOX, 0.75*d), (AL, d/2)]),
             ('fully reacted/oxidised', [(ALOX, 1.5*d)])]
    na, x = agal_alloy(d); cases.append(('dissolved in top 2 nm Ag (%.0f at%%)' % x, [(na, 2.0)]))
    for name, L in cases:
        lay = L                                     # the alloy case adds a 2 nm Ag-Al layer on top of the Ag
        e, _ = eqe(lay)
        print('%-34s %6.1f %6.2f%% %7.3f' % (name, d, 100*A_refl(lay), e))
