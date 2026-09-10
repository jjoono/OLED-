import numpy as np, json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ce_fit as cf, ce_fit2 as f2, ellipsometry_fit as ef

R = json.load(open('ce_final_result.json'))
LAB = {'#1':'VendorA ITO-1','#2':'VendorA ITO-2','#3':'VendorB IZO-1','#4':'VendorB IZO-2','#5':'VendorA ITO 2%O2'}
LOW, HIGH = 340., 1080.
fig = plt.figure(figsize=(15, 8.6))
gs = fig.add_gridspec(2, 3, hspace=0.30, wspace=0.26)

# --- (a,b) fit quality for #1 over the valid window, plus what lies outside
for col, sh in enumerate(['#1', '#3']):
    p = np.array(R[sh]['p'])
    wl, Pm, Dm = cf.load(sh)
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    Nf = f2.film_N(p, wl); Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf)); amb = np.ones_like(Nf)
    ax = fig.add_subplot(gs[0, col]); ax2 = ax.twinx()
    for a, an in enumerate(cf.ANG0 + p[2]):
        rp, rs = ef._tmm(wl, [amb, Nr, Nf, ox, si], [p[1], p[0], cf.D_OX], an)
        rr = rp/rs; ps = np.rad2deg(np.arctan(np.abs(rr))); dl = np.rad2deg(np.angle(rr)) % 360
        ax.plot(wl, Pm[:, a], color='tab:red', lw=1.1, alpha=.8)
        ax.plot(wl, ps, 'k--', lw=1.0)
        ax2.plot(wl, Dm[:, a] % 360, color='tab:green', lw=1.1, alpha=.8)
        ax2.plot(wl, dl, 'k:', lw=1.0)
    ax.axvspan(wl.min(), LOW, color='0.85'); ax.axvspan(HIGH, wl.max(), color='0.85')
    ax.set_xlim(wl.min(), wl.max()); ax.set_ylim(0, 60); ax2.set_ylim(0, 360)
    ax.set_xlabel('wavelength (nm)'); ax.set_ylabel('$\Psi$ (deg)', color='tab:red')
    ax2.set_ylabel('$\Delta$ (deg)', color='tab:green')
    ax.set_title('%s   MSE = %.1f over 340-1080 nm\n(grey = excluded: UV unmodellable / Si backside)'
                 % (LAB[sh], R[sh]['mse']), fontsize=10)

ax = fig.add_subplot(gs[0, 2]); ax.axis('off')
ax.text(0, 1, 'why the grey bands are excluded', fontsize=11, weight='bold', va='top')
ax.text(0, .88,
 '< 340 nm   a free (n,k) at every wavelength - the best any\n'
 '   homogeneous layer can do - still leaves residual 0.06-0.09\n'
 '   vs a 0.003-0.008 floor in the visible, and that is unchanged\n'
 '   for roughness 0, 3, 6, 10, 15 nm.  Structural, not fittable.\n\n'
 '> 1080 nm  Si becomes transparent at 1107 nm and the wafer\n'
 '   BACKSIDE reflection enters.  The front-surface model breaks:\n'
 '   chain k pins to 0 while n runs 1.79 -> 2.85.\n\n'
 'Fitting the full 192-1688 nm range is what produced the runaway:\n'
 'the model chases the contaminated NIR with a doubled Drude and\n'
 'then a NEGATIVE Gaussian to cancel it, giving k < 0 at 430-880 nm.',
 fontsize=8.6, va='top', family='monospace')

# --- n and k for all five
axn = fig.add_subplot(gs[1, 0]); axk = fig.add_subplot(gs[1, 1])
wlp = np.linspace(340, 1080, 500)
rows = []
for sh in ['#1', '#2', '#3', '#4', '#5']:
    N = f2.film_N(np.array(R[sh]['p']), wlp)
    axn.plot(wlp, N.real, lw=1.6, label=LAB[sh])
    axk.plot(wlp, N.imag, lw=1.6, label=LAB[sh])
    rows.append((sh, N))
axn.set_xlabel('wavelength (nm)'); axn.set_ylabel('n'); axn.grid(alpha=.3); axn.legend(fontsize=8)
axk.set_xlabel('wavelength (nm)'); axk.set_ylabel('k'); axk.grid(alpha=.3)
axk.set_ylim(bottom=0); axn.set_title('n (Gen-Osc, KK-consistent)', fontsize=10)
axk.set_title('k  (>= 0 by construction)', fontsize=10)

ax = fig.add_subplot(gs[1, 2]); ax.axis('off')
tab = 'sample            d(nm) rough  dth    n550   k550   MSE\n' + '-'*54 + '\n'
for sh, N in rows:
    p = R[sh]['p']; i = np.argmin(abs(wlp-550))
    tab += '%-16s %5.1f %5.2f %+5.2f  %5.3f  %5.3f  %5.1f\n' % (
        LAB[sh], p[0], p[1], p[2], N.real[i], N.imag[i], R[sh]['mse'])
ax.text(0, 1, tab, fontsize=9, va='top', family='monospace')
ax.text(0, .55, 'geometry is LOCKED to the transmittance-validated\nvalues; letting it float buys MSE 22 -> 16 only by\nsliding onto the n*d alias (d=46 nm, k550=0.115)\nthat the measured transmittance already rejected.',
        fontsize=8.6, va='top')
fig.suptitle('ITO / IZO 260813 - Gen-Osc model in CompleteEASE\'s frame, valid window 340-1080 nm', fontsize=12)
fig.savefig('ce_v8_check.png', dpi=125, bbox_inches='tight')
print('saved ce_v8_check.png')
for sh, N in rows:
    i = np.argmin(abs(wlp-550)); j = np.argmin(abs(wlp-633))
    print('  %-16s n550=%.3f k550=%.4f  n633=%.3f k633=%.4f' % (LAB[sh], N.real[i], N.imag[i], N.real[j], N.imag[j]))
