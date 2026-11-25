"""A physical band-edge tail for the dense perovskite: the 'zeroing and stitching' of Dasgupta et al.
(arXiv:2601.11793) with the device's own luminescence as the sub-gap probe.

Their point: ellipsometry cannot give the band-edge shape or the sub-gap extinction of a perovskite film
(roughness and scattering leave an unphysical k below the gap; forcing k to zero at a chosen wavelength
gives an arbitrary edge).  Their cure: take the band edge from a measurement that sees absorption only
(PDS or FTPS), scale it to the ellipsometric k at a stitch point just above the edge, extend the tail
log-linearly below the noise floor, and recalculate n by a (singly subtractive) Kramers-Kronig transform
anchored at the ellipsometric n of the stitch point.

No PDS exists for this film.  What does exist is the electroluminescence spectrum, and emission and
absorption of the same states are tied by detailed balance (van Roosbroeck-Shockley):

    internal emission rate(E)  ~  n^2 E^2 alpha(E) exp(-E / kT)      ->      alpha(E) ~ phi_int(E) exp(E / kT) / (n^2 E^2)

phi_int is the internal spectrum of gp_device_model (the measured EL divided by the outcoupling of the
normal direction).  Here the fit of gp_perov has k = 0 above ~565 nm by construction (no basis function
below 2.24 eV): zeroed, but with an edge whose last decade is set by the basis, not by data.

Steps
  1  alpha from the internal spectrum, 520-583 nm, in arbitrary units
  2  scaled to the ellipsometric k at the stitch wavelength (where k_SE = K_STITCH), used from there on
  3  log-linear (Urbach) extrapolation from the last stretch of usable signal
  4  n recalculated: n_SE + KK[k_stitched - k_SE].  They rebuild all of n from k and therefore need an anchor
     for the unknown constant; here only the change of k near the edge is transformed, on top of an n that
     is already Kramers-Kronig consistent, and the natural condition is that the change of n vanishes far
     from the edge - no anchor.
Above the stitch point the ellipsometric k is kept, as they do: there the luminescence-derived alpha
depends on the self-absorption correction of the device model and is not an independent measurement.

-> $ELLIPS_OUT/grahil_perov/Perov_dense_nk_stitched.csv, lib_Perov_dense_stitched.csv, stitch.json
If DASGUPTA_CSV points at their NeatPb_FaCsPbBrI_stitched_nointerp.csv, the Br = 100 % set
(Cs0.17FA0.83PbBr3) is compared: directly, and shifted along the wavelength axis so that the band edges
coincide (their own alignment step)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, csv, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_model as P, gp_perov as PV, gp_device_check as DC, gp_device_model as DM

W = _os.path.join(ELLIPS_OUT, ''); OUT = _os.path.join(ELLIPS_OUT, 'grahil_perov', '')
KT = 0.02585                      # eV, 300 K
K_STITCH = 0.10                   # stitch where the ellipsometric k has fallen to this (just above the edge)
TAIL_FIT = (566.0, 583.0)         # nm, stretch of the luminescence-derived edge used for the log-linear extrapolation
DASGUPTA = _os.environ.get('DASGUPTA_CSV', '')


def kk_delta_n(E, dk):
    """n(E) change belonging to a change dk(E) of the extinction: (2/pi) P int E' dk(E') / (E'^2 - E^2) dE'  (uniform grid, Maclaurin)"""
    h = E[1] - E[0]; out = np.zeros(len(E))
    for i in range(len(E)):
        j = np.arange((i + 1) % 2, len(E), 2)                       # every other point, skipping the pole
        out[i] = (2 / np.pi) * 2 * h * np.sum(E[j] * dk[j] / (E[j]**2 - E[i]**2))
    return out


def main():
    _os.makedirs(OUT, exist_ok=True)
    base = DC.fit_of('base'); th, Im, el = DC.measured(); el = np.interp(DM.WL, DC.WL, el)
    C = DM.components('SE', base); r = DM.evaluate(C, None, 15.0, 'lost', 'follow', (th, Im, el))
    spec = np.array(r['internal_spectrum']); wls = DM.WL
    # 1  detailed balance.  The digitised EL is power per nm: photons per eV ~ P_lambda lambda^3; alpha ~ phi_E exp(E/kT) / (n^2 E^2)
    E = P.HC / wls; n_s = PV.nk_of(base, wls).real
    alpha = spec * wls**3 * np.exp((E - E.max()) / KT) / (n_s**2 * E**2); k_rec = alpha * wls                # k = alpha lambda / 4 pi, arbitrary scale
    # fine grid
    wl = np.arange(300.0, 1000.01, 0.5); N0 = PV.nk_of(base, wl); k0 = N0.imag
    edge = (wl > 525) & (wl < 560); w_st = float(wl[edge][np.argmin(np.abs(k0[edge] - K_STITCH))])
    ok = (wls >= 520) & (wls <= TAIL_FIT[1]); lk = np.interp(wl, wls[ok], np.log(k_rec[ok]))
    scale = np.interp(w_st, wl, k0) / np.exp(np.interp(w_st, wl, lk))
    # 3  Urbach tail from the last stretch
    m = (wls >= TAIL_FIT[0]) & (wls <= TAIL_FIT[1]); pf = np.polyfit(E[m], np.log(k_rec[m] * scale), 1); e_u = 1.0 / pf[0]
    k_new = k0.copy(); a = (wl >= w_st) & (wl <= TAIL_FIT[1]); k_new[a] = scale * np.exp(lk[a])
    b = wl > TAIL_FIT[1]; k_new[b] = np.exp(np.polyval(pf, P.HC / wl[b]))
    # 4  KK on a uniform energy grid, anchored at the stitch point
    Eg = np.arange(P.HC / 1000.0, P.HC / 300.0, 0.001); wg = P.HC / Eg; dk = np.interp(wg[::-1], wl, k_new - k0)[::-1]
    dn_E = kk_delta_n(Eg, dk); dn = np.interp(wl, wg[::-1], dn_E[::-1])
    n_new = N0.real + dn
    cols = np.c_[wl, n_new, k_new, N0.real, k0]
    np.savetxt(OUT + 'Perov_dense_nk_stitched.csv', cols, delimiter=',', fmt=['%.1f', '%.5f', '%.4e', '%.5f', '%.4e'], comments='', header='wavelength_nm,n,k,n_ellipsometry,k_ellipsometry')
    np.savetxt(OUT + 'lib_Perov_dense_stitched.csv', cols[::2, :3], delimiter=',', fmt=['%.1f', '%.5f', '%.4e'], comments='', header='wavelength_nm,n,k')
    q = [520, 530, 535, 540, 545, 550, 555, 560, 570, 580, 600, 633, 700]
    tab = {str(w): dict(n=float(np.interp(w, wl, n_new)), k=float(np.interp(w, wl, k_new)), n_se=float(np.interp(w, wl, N0.real)), k_se=float(np.interp(w, wl, k0)),
                        k_luminescence=float(scale * np.exp(np.interp(w, wls[ok], np.log(k_rec[ok])))) if w <= TAIL_FIT[1] else None) for w in q}
    out = dict(stitch_wavelength_nm=w_st, k_at_stitch=K_STITCH, urbach_energy_meV=1e3 * e_u, kT_eV=KT, max_abs_dn=float(np.abs(dn).max()), table=tab,
               internal_spectrum=dict(wl=wls.tolist(), phi=spec.tolist(), peak=r['internal_peak']))
    print('stitch at %.1f nm (k = %.3f); Urbach energy %.1f meV; largest change of n %.4f' % (w_st, K_STITCH, 1e3 * e_u, np.abs(dn).max()))
    print('  nm   n_SE   k_SE     | k from luminescence (scaled) | stitched n, k')
    for w in q:
        t = tab[str(w)]; print('%5d  %.3f  %-9.3g| %-28s | %.3f  %.3g' % (w, t['n_se'], t['k_se'], ('%.3g' % t['k_luminescence']) if t['k_luminescence'] else '(extrapolated)', t['n'], t['k']))
    # ---- comparison with the published FA/Cs lead bromide set
    if DASGUPTA and _os.path.exists(DASGUPTA):
        rows = list(csv.reader(open(DASGUPTA, newline=''))); h = rows[0]; a = np.array([[float(x) if x not in ('', 'nan') else np.nan for x in rr] for rr in rows[1:]])
        wd, nd, kd = a[:, 0], a[:, h.index('n_Br_100%_FaCsPbBrI')], a[:, h.index('k_Br_100%_FaCsPbBrI')]
        half = lambda w, k: float(w[(w > 520) & (w < 580)][np.argmin(np.abs(k[(w > 520) & (w < 580)] - 0.5 * np.interp(520.0, w, k)))])
        # edge position: where k falls through 0.05 (same absolute level in both)
        lev = lambda w, k: float(w[(w > 520) & (w < 600)][np.argmin(np.abs(k[(w > 520) & (w < 600)] - 0.05))])
        shift = lev(wd, kd) - lev(wl, k_new); mt = (wd > 560) & (wd < 640) & (kd > 1e-5) & (kd < 1e-2); eu_d = 1e3 / np.polyfit(P.HC / wd[mt], np.log(kd[mt]), 1)[0]
        cmp = {str(w): dict(n_ours=float(np.interp(w, wl, n_new)), k_ours=float(np.interp(w, wl, k_new)), n_pub=float(np.interp(w, wd, nd)), k_pub=float(np.interp(w, wd, kd)),
                            n_pub_shifted=float(np.interp(w + shift, wd, nd)), k_pub_shifted=float(np.interp(w + shift, wd, kd))) for w in (350, 400, 450, 500, 520, 530, 540, 550, 560, 600, 633, 700, 800, 1000)}
        out['published'] = dict(source='Dasgupta et al., arXiv:2601.11793, Cs0.17FA0.83PbBr3 (Br 100 %), stitched', edge_shift_nm=shift, urbach_energy_meV=float(eu_d), table=cmp)
        print('published Cs0.17FA0.83PbBr3: band edge %.1f nm to the red of this film; its Urbach energy %.1f meV' % (shift, eu_d))
        print('  nm   ours n, k        | published n, k     | published shifted by the edge difference')
        for w, t in cmp.items(): print('%5s  %.3f %-9.3g | %.3f %-9.3g | %.3f %-9.3g' % (w, t['n_ours'], t['k_ours'], t['n_pub'], t['k_pub'], t['n_pub_shifted'], t['k_pub_shifted']))
    json.dump(out, open(OUT + 'stitch.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
