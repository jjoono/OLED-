"""Parameter table, n,k CSVs and a figure for the dispersions that went into the .mod files."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, csv, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_genosc as G, tl_make_mod as MM
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

W = _os.path.join(ELLIPS_OUT, '')
OUT = _os.path.join(ELLIPS_OUT, 'tooling_mod', '')
HC = 1239.841984

def main():
    CE = json.load(open(W + 'ce_genosc.json'))
    use = {'NPB': 'pick', 'TAPC': 'pick', 'Bphen': 'with_roughness'}
    wl = np.arange(250.0, 1001.0, 1.0)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6)); col = {'NPB': '#1f77b4', 'TAPC': '#d62728', 'Bphen': '#2ca02c'}
    for name, key in use.items():
        rec = CE[name][key]; p = np.array(rec['p']); ci = np.array(rec['ci']); free = MM.lock_flags(rec)
        print('%s  (Einf 1 fixed, %d Gaussians, MSE %s)' % (name, (len(p) - 2) // 3, ', '.join('%.2f' % m for m in rec['mse_per'])))
        print('   UV Pole  Amp %7.2f +-%5.2f eV^2   En %5.2f +-%4.2f eV   %s' % (p[0], ci[0], p[1], ci[1], '' if free[0] and free[1] else '(locked: ' + ', '.join(n for n, f in (('Amp', free[0]), ('En', free[1])) if not f) + ')'))
        for g in range((len(p) - 2) // 3):
            sl = slice(2 + 3 * g, 5 + 3 * g); a, b, e = p[sl]; ca, cb, ce = ci[sl]
            lk = [n for n, f in zip(('Amp', 'Br', 'En'), free[sl]) if not f]
            print('   G%d  Amp %5.3f +-%5.3f   Br %5.3f +-%5.3f eV   En %5.3f +-%5.3f eV (%3.0f nm)  %s'
                  % (g + 1, a, ca, b, cb, e, ce, HC / e, '(locked: ' + ', '.join(lk) + ')' if lk else ''))
        N = G.N_of(p, wl)
        with open(OUT + 'nk_%s_GenOsc.csv' % name, 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['wavelength_nm', 'n', 'k'])
            for a, z in zip(wl, N): w.writerow(['%.0f' % a, '%.5f' % z.real, '%.5f' % max(z.imag, 0.0)])
        print('   n,k: ' + '  '.join('%3.0f: %.3f+%.4fj' % (x, N[int(x - 250)].real, N[int(x - 250)].imag) for x in (350, 400, 450, 550, 633)))
        ax[0].plot(wl, N.real, color=col[name], label=name); ax[1].plot(wl, N.imag, color=col[name], label=name)
    lib = np.genfromtxt(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', 'dft_seedlayer_screen', 'data', 'nk', 'NPB.csv'), delimiter=',', names=True)
    ax[0].plot(lib['wavelength_nm'][::15], lib['n'][::15], 'o', ms=3, mfc='none', color=col['NPB'], label='lab NPB.csv')
    for a, lab in zip(ax, ('n', 'k')):
        a.set_xlabel('wavelength (nm)'); a.set_ylabel(lab); a.grid(alpha=.3); a.legend(fontsize=8); a.set_xlim(250, 800)
    ax[0].set_title('Gen-Osc (Einf 1 + UV pole + Gaussians), CompleteEASE frame'); ax[1].set_title('k: zero below each absorption onset')
    fig.tight_layout(); fig.savefig(OUT + 'nk_GenOsc.png', dpi=140)

if __name__ == '__main__':
    main()
