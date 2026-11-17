import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, csv, os, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_osc as O
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = _os.path.join(ELLIPS_OUT, '')
OUT = _os.path.join(ELLIPS_OUT, 'tooling_nk', ''); _os.makedirs(OUT, exist_ok=True)
OSC = json.load(open(W + 'osc.json')); SEP = json.load(open(W + 'osc_npb_sep.json'))
SI = json.load(open(W + 'si_id.json')); GG = json.load(open(W + 'glass_global.json'))
GUESS = {'1_21': ('NPB', 30), '1_31': ('NPB', 30), '1_22': ('NPB', 50), '1_32': ('NPB', 50),
         '1_23': ('NPB', 10), '1_33': ('NPB', 10), '1_24': ('TAPC', 50), '1_34': ('TAPC', 50),
         '2_21': ('Bphen', 50), '2_22': ('Bphen', 50), '2_31': ('Bphen', 50), '2_32': ('Bphen', 50),
         '2_23': ('Bphen', 10), '2_24': ('Bphen', 10), '2_33': ('Bphen', 10), '2_34': ('Bphen', 10)}
CONT = {'2_23', '2_24', '2_33', '2_34'}     # not a continuous film

# ---------------- n,k tables
wl = np.arange(250.0, 1001.0, 1.0)
disp = {'NPB': OSC['NPB']['q'], 'TAPC': OSC['TAPC']['q'], 'Bphen': OSC['Bphen']['q'],
        'NPB_50nm_own': SEP['NPB_22']['q'], 'NPB_30nm_own': SEP['NPB_21']['q']}
for name, q in disp.items():
    N = O.N_of(np.array(q), wl)
    with open(OUT + 'nk_%s.csv' % name, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['wavelength_nm', 'n', 'k'])
        for a, z in zip(wl, N): w.writerow(['%.0f' % a, '%.5f' % z.real, '%.5f' % max(z.imag, 0.0)])
g = GG['normalised']['x']
with open(OUT + 'n_glass_substrate.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['wavelength_nm', 'n'])
    for a in wl: w.writerow(['%.0f' % a, '%.5f' % (g[0] + g[1] / (a / 1000)**2)])

# ---------------- thickness table
SHG = ['1_31', '1_32', '1_33', '1_34', '2_31', '2_32', '2_33', '2_34']
rows = []
for sh in sorted(GUESS):
    mat_g, nom = GUESS[sh]
    if sh.split('_')[1].startswith('2'):
        t = SI[sh]; b = t['best']; d = t[b][0]; mse = t[b][1]
        lo, hi = sorted([t['d_ox1'], t['d_ox3']]); sub = 'Si'; how = 'native oxide 1-3 nm'
        others = sorted(t[m][1] for m in ('NPB', 'TAPC', 'Bphen') if m != b)[0]
    else:
        j = SHG.index(sh); b = mat_g; d = GG['normalised']['x'][3 + j]; mse = GG['normalised']['per'][j]
        ds = [GG[c]['x'][3 + j] for c in GG]; lo, hi = min(ds), max(ds); sub = 'glass'; how = 'depolarisation convention'
        others = np.nan
    note = ''
    if sh in CONT:
        b = 'undetermined'
        if sub == 'Si':
            note = 'no homogeneous film fits (best MSE ~39 for any material, any oxide, angle offset or void fraction); thickness not reported'
            d = lo = hi = np.nan
        else:
            note = 'near-bare glass; value is the Bphen-equivalent thickness of whatever is there, not a film thickness'
    rows.append(dict(sheet=sh, substrate=sub, guess=mat_g, nominal_nm=nom, identified=b,
                     d_nm=d, d_lo=lo, d_hi=hi, ratio=d / nom, mse=mse, next_material_mse=others,
                     spread_from=how, note=note))
with open(OUT + 'thickness_summary.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: ('' if isinstance(v, float) and np.isnan(v) else '%.3f' % v if isinstance(v, float) else v) for k, v in r.items()})
for r in rows:
    print('%-5s %-5s %-6s %3d | %-24s %6.2f nm [%5.2f-%5.2f]  ratio %.3f  MSE %5.1f  next %6.1f'
          % (r['sheet'], r['substrate'], r['guess'], r['nominal_nm'], r['identified'], r['d_nm'], r['d_lo'], r['d_hi'],
             r['ratio'], r['mse'], r['next_material_mse']))

# ---------------- figures
lib = np.genfromtxt(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..', 'dft_seedlayer_screen', 'data', 'nk', 'NPB.csv'), delimiter=',', names=True)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
col = {'NPB': '#1f77b4', 'TAPC': '#d62728', 'Bphen': '#2ca02c'}
for m in ('NPB', 'TAPC', 'Bphen'):
    N = O.N_of(np.array(disp[m]), wl)
    ax[0].plot(wl, N.real, color=col[m], label=m); ax[1].plot(wl, N.imag, color=col[m], label=m)
N = O.N_of(np.array(disp['NPB_30nm_own']), wl)
ax[0].plot(wl, N.real, color=col['NPB'], ls=':', lw=1, label='NPB, 30 nm piece alone')
ax[1].plot(wl, N.imag, color=col['NPB'], ls=':', lw=1)
ax[0].plot(lib['wavelength_nm'][::15], lib['n'][::15], 'o', ms=3, color=col['NPB'], mfc='none', label='lab library NPB.csv')
for a, lab in zip(ax, ('n', 'k')):
    a.set_xlabel('wavelength (nm)'); a.set_ylabel(lab); a.grid(alpha=.3); a.legend(fontsize=8); a.set_xlim(250, 800)
ax[0].set_title('refractive index (fitted on the Si pieces)'); ax[1].set_title('extinction coefficient')
fig.tight_layout(); fig.savefig(OUT + 'nk_NPB_TAPC_Bphen.png', dpi=140); plt.close(fig)

fig, ax = plt.subplots(figsize=(7.5, 4.6))
pairs = [('1_21', '1_31'), ('1_22', '1_32'), ('1_23', '1_33'), ('1_24', '1_34'), ('2_21', '2_31'), ('2_22', '2_32')]
x = np.arange(len(pairs))
dsi = [next(r for r in rows if r['sheet'] == a)['d_nm'] for a, b in pairs]
dgl = [next(r for r in rows if r['sheet'] == b)['d_nm'] for a, b in pairs]
nom = [GUESS[a][1] for a, b in pairs]
ax.bar(x - 0.2, dsi, 0.38, label='Si piece', color='#555'); ax.bar(x + 0.2, dgl, 0.38, label='glass piece', color='#9ab')
ax.scatter(x, nom, marker='_', s=900, color='k', zorder=3, label='nominal (your guess)')
ax.set_xticks(x); ax.set_xticklabels(['%s / %s\n%s' % (a, b, GUESS[a][0]) for a, b in pairs], fontsize=8)
ax.set_ylabel('fitted thickness (nm)'); ax.grid(axis='y', alpha=.3); ax.legend(fontsize=8)
ax.set_title('continuous films: Si and glass pieces vs nominal')
fig.tight_layout(); fig.savefig(OUT + 'thickness_si_vs_glass.png', dpi=140); plt.close(fig)
print('->', OUT, sorted(os.listdir(OUT)))
