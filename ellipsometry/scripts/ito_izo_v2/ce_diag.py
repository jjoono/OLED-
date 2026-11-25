import numpy as np, ce_fit as cf, ellipsometry_fit as ef

wl, Pm, Dm = cf.load('#1')
# my earlier solution, expressed in this parameter vector (no UV Gaussian)
p_old = np.array([51., 3., 0., 2.5698279, 1.0644047e-3, 6.582119569,
                  165.546637, 1.05505954, 4.0, 3.33679487,
                  0.0, 2.5, 6.0, 0.07201886, 0.82494932, 1.91114700])

def run(tag, N_ox, N_si, lo, hi, p):
    m = (wl >= lo) & (wl <= hi)
    r = cf.make_resid(wl[m], Pm[m], Dm[m], N_ox[m], N_si[m])
    res = r(p)
    print('  %-34s %4.0f-%4.0f nm  MSE = %8.2f' % (tag, lo, hi, cf.mse(res, 16)))
    return res

ce_ox, ce_si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
z = np.load('izo_data.npz')
my_si = np.sqrt((np.interp(wl, z['si_wl'], z['si_e1']) + 1j*np.interp(wl, z['si_wl'], z['si_e2'])).astype(complex))
my_ox = np.sqrt((np.interp(wl, z['ox_wl'], z['ox_e1']) + 1j*np.interp(wl, z['ox_wl'], z['ox_e2'])).astype(complex))

print('my old solution, evaluated with each substrate set:')
run('CompleteEASE SI_JAW+NTVE_JAW', ce_ox, ce_si, 251, 1078, p_old)
run('my own Si + SiO2',             my_ox, my_si, 251, 1078, p_old)
run('CompleteEASE SI_JAW+NTVE_JAW', ce_ox, ce_si, 251, 1688, p_old)
run('CompleteEASE SI_JAW+NTVE_JAW', ce_ox, ce_si, 400,  900, p_old)

print('\nnative-oxide index at 633 nm:  NTVE_JAW n=%.4f   mine n=%.4f'
      % (np.interp(633, wl, ce_ox.real), np.interp(633, wl, my_ox.real)))
print('substrate  n,k at 633 nm:      SI_JAW %.4f%+.4fj   mine %.4f%+.4fj'
      % (np.interp(633, wl, ce_si.real), np.interp(633, wl, ce_si.imag),
         np.interp(633, wl, my_si.real), np.interp(633, wl, my_si.imag)))

print('\nmeasured vs my-old-model Psi/Delta (CompleteEASE substrates, 65 deg):')
C = cf.forward(p_old, wl, ce_ox, ce_si)
n0, c0, s0 = C[0]
psi_c = np.rad2deg(0.5*np.arccos(np.clip(n0, -1, 1)))
del_c = np.rad2deg(np.arctan2(s0, c0))
for t in [250, 400, 600, 800, 1000, 1200, 1500, 1688]:
    i = np.argmin(abs(wl-t))
    print('  %6.1f nm   Psi meas %7.3f  calc %7.3f    Delta meas %8.3f  calc %8.3f'
          % (wl[i], Pm[i,0], psi_c[i], Dm[i,0], del_c[i]))
