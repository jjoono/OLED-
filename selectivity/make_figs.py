"""Figures + Excel source data for the directional-selectivity study (figSel1-4)."""
import os, json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
import sel, common, bound

os.makedirs('figures', exist_ok=True)
plt.rcParams.update({'font.size': 7, 'axes.titlesize': 7.5, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.color': '#e4e3df', 'grid.linewidth': 0.5, 'lines.linewidth': 1.4})
SW = common.Sweep()
DV = dict(np.load('data/devices.npz')); SL = dict(np.load('data/slabs.npz'))
MLA = dict(np.load('data/ref_mla.npz')); CAV = np.load('data/ref_cavity.npz')
RND = json.load(open('data/random.json')); KITO = json.load(open('data/kito.json'))
OPT = json.load(open('data/optfilter.json')) if os.path.exists('data/optfilter.json') else []
ORDER = ['Al', 'Ag', 'DBRrear', 'LL15', 'LL18']
COL = dict(zip(ORDER, ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']))
MRK = dict(zip(ORDER, ['o', 's', 'D', '^', 'v']))
LAB = {'Al': 'Al, glass', 'Ag': 'Ag, glass', 'DBRrear': 'rear DBR, glass', 'LL15': 'low-loss, n 1.51', 'LL18': 'low-loss, n 1.80'}
INK2 = '#52514e'
G30 = common.SPECTRA['G30']


def Aprime(st, spec=G30):
    return float(spec @ (1 - (DV[f'{st}_R'] * sel.WL).sum(1)))


def esub(st, spec=G30):
    return float(spec @ DV[f'{st}_eta_sub'])


def fwd(e):
    return e[common.BI[5.0]] / np.sin(np.radians(5.0)) ** 2


def save(fig, name, sheets):
    fig.savefig(f'figures/{name}.png', dpi=300, bbox_inches='tight'); fig.savefig(f'figures/{name}.pdf', bbox_inches='tight')
    wb = Workbook(); wb.remove(wb.active)
    for title, header, rows in sheets:
        ws = wb.create_sheet(title[:31]); ws.append(header)
        for r in rows:
            ws.append([float(x) if isinstance(x, (np.floating, np.integer)) else x for x in r])
    wb.save(f'figures/{name}.xlsx'); plt.close(fig); print('saved', name)


def layers_from(ds, mat):
    H, L = sel.DBR_MATS[mat]
    return ('layers', [(H if k % 2 == 0 else L, d) for k, d in enumerate(ds) if d > 0.05])


def eval_filter(st, spec, filt, rand='full'):
    dev = dict(eta_sub=DV[f'{st}_eta_sub'], P=DV[f'{st}_P'], R=DV[f'{st}_R']); ns = float(DV[f'{st}_n_sub'])
    Rf, Tf = sel.filter_RT(filt, ns)
    Ts, Rs = sel.mixer(1.0) if rand == 'full' else tuple(SL['g0.8_t1.0_w1'])
    return sel.metrics(dev['eta_sub'], sel.cones(sel.solve(dev, Rf, Tf, Ts, Rs), ns), spec)


# =========================================================== Fig Sel1: ceiling vs round-trip loss
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw=dict(wspace=0.42))
A = np.logspace(np.log10(0.005), np.log10(0.3), 200); sh = []
BC = {10: '#2a78d6', 20: '#eb6834', 30: '#1baf7a', 45: '#6250d6'}
for b, c in BC.items():
    for n, ls in ((1.51, '-'), (1.8, '--')):
        y = bound.eta_cone(A, b, n)
        ax[0].plot(A, y, ls, color=c, lw=1.3, label=f'{b}° cone' if n == 1.51 else None)
        sh.append((f'a_ceiling_b{b}_n{n}', ['A_prime', 'eta_cone_max'], list(zip(A, y))))
pts = []
for st in ORDER:
    for b in (10, 20, 30):
        e = SW.ref(st, 'full', G30, ('ideal', float(b)))[common.BI[float(b)]] / esub(st)
        ax[0].plot(Aprime(st), e, MRK[st], ms=4, mfc='white', mec=COL[st], mew=1.1)
        pts.append([st, Aprime(st), b, e])
for st in ORDER:
    ax[0].plot([], [], MRK[st], mfc='white', mec=COL[st], label=LAB[st])
sh.append(('a_devices_ideal_filter', ['stack', 'A_prime_G30', 'beta', 'eta_cone_sta (matrix model)'], pts))
ax[0].set_xscale('log'); ax[0].set_xlabel("round-trip loss A'"); ax[0].set_ylabel('fraction of substrate light into cone')
ax[0].set_title('(a) ceiling, ideal filter + randomizer', loc='left'); ax[0].legend(fontsize=4.8, ncol=3, loc='upper left', bbox_to_anchor=(0.0, -0.22), frameon=False)
ax[0].text(0.31, 0.93, 'solid n 1.51, dashed n 1.80', fontsize=5, color=INK2, transform=ax[0].transAxes)
for b, c in list(BC.items())[:3]:
    y = bound.fwd_gain(A, b, 1.51)
    ax[1].plot(A, y, color=c, label=f'{b}° cone'); sh.append((f'b_fwd_gain_b{b}', ['A_prime', 'forward_luminance_gain'], list(zip(A, y))))
for st in ORDER:
    ax[1].axvline(Aprime(st), color=COL[st], lw=0.6, ls=':')
ax[1].set_xscale('log'); ax[1].set_yscale('log'); ax[1].set_xlabel("round-trip loss A'"); ax[1].set_ylabel('on-axis luminance gain vs no filter')
ax[1].set_title('(b) forward gain ceiling (n 1.51)', loc='left'); ax[1].legend(fontsize=5, frameon=False)
k = [r['k'] for r in KITO]; o = [r['G30'] for r in KITO]
nof = np.array([x['nofilter']['EQE_total'] for x in o]); i20 = np.array([x['ideal20']['EQE_20'] for x in o]); d20 = np.array([x['dbr20']['EQE_20'] for x in o])
for y, c, lab in ((nof, '#7d7c76', 'total EQE, no filter'), (i20, '#eb6834', '20° cone, ideal filter'), (d20, '#2a78d6', '20° cone, edge DBR')):
    ax[2].plot(k, y / y[0], 'o-', color=c, ms=3, label=lab)
ax2 = ax[2].secondary_xaxis('top', functions=(lambda x: np.interp(x, k, [r['G30']['A_prime'] for r in KITO]),
                                               lambda a: np.interp(a, [r['G30']['A_prime'] for r in KITO], k)))
ax2.set_xlabel("A'", fontsize=6); ax2.tick_params(labelsize=5.5); ax2.set_xticks([0.02, 0.03, 0.04, 0.05, 0.06])
ax[2].set_xlabel('ITO extinction coefficient k (550 nm)'); ax[2].set_ylabel('relative to k = 0')
ax[2].set_title('(c) low-loss stack, n 1.80', loc='left', pad=14); ax[2].legend(fontsize=5, frameon=False, loc='lower left')
sh.append(('c_kITO', ['k_ITO_550', 'A_prime', 'EQE_total_nofilter', 'EQE20_ideal', 'EQE20_edgeDBR'],
           [[r['k'], r['G30']['A_prime'], a, b, c] for r, a, b, c in zip(KITO, nof, i20, d20)]))
save(fig, 'figSel1_ceiling', sh)

# =========================================================== Fig Sel2: realistic filters vs bandwidth
FW = [10, 20, 30, 45, 60, 80, 100]; REAL = ['blue QD', 'v-DABNA (HF)', 'Ir(ppy)2acac', 'Ir(dmppy-ph)2tmd']
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.42)); sh = []
for k_, st in enumerate(('LL18', 'Al')):
    a = ax[k_]; rows = []
    ide = [SW.ref(st, 'full', common.SPECTRA[f'G{f}'], ('ideal', 20.0))[common.BI[20.0]] for f in FW]
    nf = [SW.ref(st, 'full', common.SPECTRA[f'G{f}'], ('fresnel',))[common.BI[20.0]] for f in FW]
    ml = [common.SPECTRA[f'G{f}'] @ (DV[f'{st}_eta_sub'][:, None] * MLA[st])[:, common.BI[20.0]] for f in FW]
    a.plot(FW, ide, color='#0b0b0b', lw=1.0, ls='--', label='ideal filter')
    a.plot(FW, nf, color='#7d7c76', lw=1.0, label='scatterer, no filter')
    a.plot(FW, ml, color='#b5b4ad', lw=1.0, ls='-.', label='hemispherical MLA')
    for (mat, N), c in {('TiO2/SiO2', 8): '#8fb8ea', ('TiO2/SiO2', 16): '#2a78d6', ('ZnS/LiF', 16): '#1baf7a'}.items():
        y = [SW.best_edge(st, 'full', common.SPECTRA[f'G{f}'], 20.0, mat=mat, N=N)[1][common.BI[20.0]] for f in FW]
        a.plot(FW, y, 'o-', color=c, ms=2.5, label=f'edge DBR {mat}, {N} pairs'); rows.append((f'{mat} N{N}', y))
    for sp, mk in zip(REAL, ('P', 'X', 'p', 'h')):
        e = SW.best_edge(st, 'full', common.SPECTRA[sp], 20.0, mat='TiO2/SiO2', N=16)[1][common.BI[20.0]]
        a.plot(common.FWHM[sp], e, mk, color='#2a78d6', ms=4.5, mec='#0b0b0b', mew=0.3)
        off = {'blue QD': (4, -8), 'v-DABNA (HF)': (4, 3), 'Ir(ppy)2acac': (-30, 6), 'Ir(dmppy-ph)2tmd': (3, -9)}[sp]
        a.annotate(sp.split(' (')[0], (common.FWHM[sp], e), xytext=off, textcoords='offset points', fontsize=4.6, color=INK2)
        rows.append((sp + ' (TiO2/SiO2 16, at its FWHM)', [e]))
    for o in OPT:
        if o['stack'] == st and o['beta'] == 20.0 and o['spec'] in common.SPECTRA:
            e = eval_filter(st, common.SPECTRA[o['spec']], layers_from(o['ds'], o['mat']))['EQE_20']
            a.plot(common.FWHM[o['spec']], e, '*', color='#eb6834', ms=8, mec='#0b0b0b', mew=0.4, zorder=5)
            rows.append((f"optimized {o['mat']} {o['M']} layers, {o['spec']}", [e]))
    a.plot([], [], '*', color='#eb6834', ms=7, mec='#0b0b0b', mew=0.4, label='layer-optimized filter')
    a.set_xlabel('emission FWHM (nm, Gaussian at 530 nm)'); a.set_ylabel('EQE into 20° cone')
    a.set_title(f"({'ab'[k_]}) {LAB[st]} (A' = {Aprime(st):.3f})", loc='left'); a.set_ylim(0, None)
    if k_ == 0:
        a.legend(fontsize=4.6, frameon=False, loc='upper right', bbox_to_anchor=(1.0, 0.88))
    sh.append((f"{'ab'[k_]}_{st}", ['curve'] + [f'FWHM {f}' for f in FW],
               [['ideal'] + ide, ['no filter'] + nf, ['MLA'] + ml] + [[n] + list(v) for n, v in rows]))
rows = []
for st in ORDER:
    fr = [SW.best_edge(st, 'full', common.SPECTRA[f'G{f}'], 20.0, mat='TiO2/SiO2')[1][common.BI[20.0]] /
          SW.ref(st, 'full', common.SPECTRA[f'G{f}'], ('ideal', 20.0))[common.BI[20.0]] for f in FW]
    ax[2].plot(FW, fr, MRK[st] + '-', color=COL[st], ms=3, label=f"{LAB[st]} (A' {Aprime(st):.3f})"); rows.append([st, Aprime(st)] + fr)
ax[2].set_xlabel('emission FWHM (nm)'); ax[2].set_ylabel('edge DBR / ideal filter (20° cone)'); ax[2].set_ylim(0, 1)
ax[2].set_title('(c) fraction of the ceiling realized', loc='left'); ax[2].legend(fontsize=4.8, frameon=False, loc='lower left')
sh.append(('c_realized_fraction', ['stack', 'A_prime'] + [f'FWHM {f}' for f in FW], rows))
save(fig, 'figSel2_real_filters', sh)

# =========================================================== Fig Sel3: randomization
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), gridspec_kw=dict(wspace=0.42)); sh = []
get = lambda st, kind, sp='G30': next(r for r in RND if r['stack'] == st and r['spec'] == sp and r['design'].startswith(kind))
MS = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 1.0]
for st in ('LL18', 'Al'):
    for kind, ls, lab in (("('ideal'", '--', 'ideal filter'), ("('edge'", '-', 'edge DBR'), ("('fresnel'", ':', 'no filter')):
        r = get(st, kind)['res']; y = [r[f'm{m}']['EQE_20'] for m in MS]
        ax[0].plot(MS, y, ls, color=COL[st], marker=MRK[st], ms=2.5, label=f'{LAB[st]}, {lab}')
        sh.append((f'a_{st}_{lab.replace(" ", "_")}', ['m', 'EQE_20'], list(zip(MS, y))))
ax[0].set_xlabel('randomized fraction per traversal, m'); ax[0].set_ylabel('EQE into 20° cone')
ax[0].set_title('(a) partial randomization', loc='left'); ax[0].legend(fontsize=4.4, frameon=False, ncol=1)
TAU = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0]; GC = {0.0: '#2a78d6', 0.5: '#1baf7a', 0.8: '#eb6834', 0.9: '#6250d6'}
r = get('LL18', "('ideal'")['res']; full = r['m1.0']['EQE_20']
for g, c in GC.items():
    y = [r[f'g{g}_t{t}_w1']['EQE_20'] / full for t in TAU]
    ax[1].plot(TAU, y, 'o-', color=c, ms=2.5, label=f'g = {g}'); sh.append((f'b_g{g}', ['tau_transport', 'EQE20/full_randomizer'], list(zip(TAU, y))))
ax[1].set_xscale('log'); ax[1].set_xlabel('transport optical thickness τ* = S(1−g)'); ax[1].set_ylabel('relative to full randomization')
ax[1].set_title('(b) scattering slab, low-loss n 1.80, ideal filter', loc='left', fontsize=6.5); ax[1].legend(fontsize=5, frameon=False); ax[1].set_ylim(0, 1.05)
W = ['w1', 'w0.9995', 'w0.999', 'w0.995']; WX = [0, 5e-4, 1e-3, 5e-3]; rows = []
for st in ORDER:
    rr = get(st, "('ideal'")['res']; y = [rr[f'g0.8_t1.0_{w}']['EQE_20'] / rr['g0.8_t1.0_w1']['EQE_20'] for w in W]
    ax[2].plot(WX, y, MRK[st] + '-', color=COL[st], ms=3, label=LAB[st]); rows.append([st, Aprime(st)] + y)
    rf = get(st, "('fresnel'")['res']; yf = [rf[f'g0.8_t1.0_{w}']['EQE_20'] / rf['g0.8_t1.0_w1']['EQE_20'] for w in W]
    ax[2].plot(WX, yf, MRK[st] + ':', color=COL[st], ms=2, mfc='white', lw=0.8); rows.append([st + ' (no filter)', Aprime(st)] + yf)
ax[2].set_xlabel('scatterer absorption per event, 1 − albedo'); ax[2].set_ylabel('EQE (20° cone) relative to lossless')
ax[2].set_title('(c) absorption in the randomizer', loc='left'); ax[2].legend(fontsize=4.8, frameon=False, loc='lower left')
ax[2].text(0.97, 0.97, 'solid: ideal filter\ndotted: no filter', ha='right', va='top', fontsize=5, color=INK2, transform=ax[2].transAxes)
sh.append(('c_albedo', ['stack', 'A_prime', 'w=1', 'w=0.9995', 'w=0.999', 'w=0.995'], rows))
save(fig, 'figSel3_randomization', sh)

# =========================================================== Fig Sel4: comparison
CATS = ['planar', 'MLA', 'strong cavity', 'scatterer only', 'scatterer + edge DBR', 'scatterer + optimized filter', 'ideal filter']
CC = ['#b5b4ad', '#7d7c76', '#6250d6', '#1baf7a', '#2a78d6', '#eb6834', '#ffffff']


def cav_best(spec, beta, cath):
    E = np.einsum('l,kl,klb->kb', spec, CAV['eta_sub'], CAV['cones'])
    ks = [i for i, k in enumerate(CAV['keys']) if k.startswith(cath + '|')]
    col = common.BI[5.0] if beta == 5.0 else common.BI[beta]
    return E[max(ks, key=lambda i: E[i, col])]


fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.25)); rows = []
for p_, metric in enumerate(('EQE into 20° cone', 'on-axis luminance (Lambertian-equivalent EQE)')):
    for i, st in enumerate(ORDER):
        b = 20.0
        v = {'planar': SW.ref(st, 'none', G30, ('fresnel',)), 'MLA': G30 @ (DV[f'{st}_eta_sub'][:, None] * MLA[st]),
             'scatterer only': SW.ref(st, 'scat', G30, ('fresnel',)),
             'scatterer + edge DBR': SW.best_edge(st, 'scat', G30, b, key='cone' if p_ == 0 else 'fwd')[1],
             'ideal filter': SW.ref(st, 'scat', G30, ('ideal', b if p_ == 0 else 5.0))}
        if st in ('Al', 'Ag'):
            v['strong cavity'] = cav_best(G30, b if p_ == 0 else 5.0, st)
        oo = [o for o in OPT if o['stack'] == st and o['spec'] == 'G30' and o['beta'] == 20.0 and o['mat'] == 'TiO2/SiO2']
        if oo:
            m = eval_filter(st, G30, layers_from(oo[0]['ds'], oo[0]['mat']), rand='scat')
            v['scatterer + optimized filter'] = np.array([m[f'EQE_{bb:g}'] for bb in sel.BETAS])
        for j, cat in enumerate(CATS):
            if cat not in v:
                continue
            y = v[cat][common.BI[b]] if p_ == 0 else fwd(v[cat])
            x = i + (j - 3) * 0.12
            ax[p_].bar(x, y, 0.11, color=CC[j], edgecolor='#0b0b0b' if cat == 'ideal filter' else 'none', lw=0.6,
                       hatch='////' if cat == 'ideal filter' else None, label=cat if i == 0 or (cat == 'strong cavity' and i == 0) else None)
            rows.append([metric, st, cat, y])
    ax[p_].set_xticks(range(len(ORDER))); ax[p_].set_xticklabels([LAB[s].replace(', ', '\n') for s in ORDER], fontsize=5.5)
    ax[p_].set_ylabel(metric); ax[p_].grid(axis='x', visible=False)
    ax[p_].set_title(['(a) power into a 20° cone', '(b) on-axis luminance (ideal: 5° cone)'][p_], loc='left')
h, l = ax[0].get_legend_handles_labels(); hs = dict(zip(l, h))
fig.legend([hs[c] for c in CATS if c in hs], [c for c in CATS if c in hs], ncol=7, fontsize=5, frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.06))
save(fig, 'figSel4_comparison', [('comparison', ['metric', 'stack', 'configuration', 'value'], rows)])

# =========================================================== optimized filters: full-model check
import csv
with open('results/optimized_filters.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['stack', 'spectrum', 'beta', 'material', 'layers', 'EQE_beta_fast', 'EQE_beta_full_randomizer', 'EQE_beta_slab_g0.8_tau1',
                'EQE_total_slab', 'EQE_fwd_slab', 'ideal_EQE_beta_slab', 'edge16_EQE_beta_slab', 'total_thickness_nm', 'thicknesses_nm'])
    for o in OPT:
        spec = common.SPECTRA[o['spec']]; flt = layers_from(o['ds'], o['mat']); b = o['beta']
        mf = eval_filter(o['stack'], spec, flt, 'full'); ms = eval_filter(o['stack'], spec, flt, 'scat')
        ide = SW.ref(o['stack'], 'scat', spec, ('ideal', b))[common.BI[b]]
        edg = SW.best_edge(o['stack'], 'scat', spec, b, mat=o['mat'], N=16)[1][common.BI[b]]
        w.writerow([o['stack'], o['spec'], b, o['mat'], o['M'], o['eqe'], mf[f'EQE_{b:g}'], ms[f'EQE_{b:g}'], ms['EQE_total'], ms['EQE_fwd'],
                    ide, edg, sum(o['ds']), ' '.join(f'{d:.1f}' for d in o['ds'])])
print(open('results/optimized_filters.csv').read())

# =========================================================== Fig Sel5: display metrics beyond on-axis efficiency
if os.path.exists('data/display_refs.json') and os.path.exists('data/display_dbrcav.json'):
    REF = json.load(open('data/display_refs.json')); DC = json.load(open('data/display_dbrcav.json')); ang = np.array(REF['angles'])
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw=dict(wspace=0.4)); sh = []
    sp = 'G30'
    pick = [('Al planar (ITO 150)', 'Al planar', '#b5b4ad'), ('low-loss n1.8 + MLA', 'low-loss + MLA', '#7d7c76'),
            ('thin-Ag strong cavity', 'thin-Ag cavity', '#6250d6'), ('DBR cavity N6', 'DBR cavity, 6 pairs', '#1baf7a'),
            ('DBR cavity N10', 'DBR cavity, 10 pairs', '#eda100'), ('scatterer + edge DBR', 'low-loss + scatterer + edge DBR', '#2a78d6')]
    allc = {**{k.split('|', 1)[1]: v for k, v in REF['configs'].items() if k.startswith(sp + '|')},
            **{k.split('|', 1)[1]: v for k, v in DC.items() if k.startswith(sp + '|')}}
    for key, lab, c in pick:
        k = next((kk for kk in allc if key in kk), None)
        if k is None:
            continue
        m = allc[k]
        ax[0].plot(ang, np.array(m['lum_rel']), color=c, label=lab); ax[1].plot(ang, m['duv'], color=c, label=lab)
        sh.append((lab[:28], ['angle_deg', 'luminance_rel', 'du_v'], list(zip(ang, m['lum_rel'], m['duv']))))
    ax[0].set_xlim(0, 80); ax[0].set_xlabel('viewing angle (deg)'); ax[0].set_ylabel('luminance / on-axis luminance')
    ax[0].set_title('(a) angular luminance, 530/30 nm', loc='left'); ax[0].legend(fontsize=4.8, frameon=False)
    ax[1].set_xlim(0, 60); ax[1].set_xlabel('viewing angle (deg)'); ax[1].set_ylabel("colour shift Δu'v' from on-axis")
    ax[1].axhline(0.02, color=INK2, lw=0.6, ls=':'); ax[1].text(1, 0.021, "0.02", fontsize=5, color=INK2)
    ax[1].set_title('(b) colour shift, 530/30 nm', loc='left')
    rows = []
    for spn, mk in (('G30', 'o'), ('Ir(ppy)2acac', 's'), ('G10', '^')):
        cc = {**{k.split('|', 1)[1]: v for k, v in REF['configs'].items() if k.startswith(spn + '|')},
              **{k.split('|', 1)[1]: v for k, v in DC.items() if k.startswith(spn + '|')}}
        for key, lab, c in pick + [('ideal 20 deg', 'ideal 20° filter', '#0b0b0b')]:
            k = next((kk for kk in cc if key in kk), None)
            if k is None:
                continue
            m = cc[k]; ax[2].plot(m['EQE'], m['cdA0'], mk, color=c, ms=4, mec='white', mew=0.4); rows.append([spn, lab, m['EQE'], m['cdA0'], m['duv30'], m['duv60'], m['half_lum_angle']])
    for spn, mk in (('G30', 'o'), ('Ir(ppy)2acac', 's'), ('G10', '^')):
        ax[2].plot([], [], mk, color='#7d7c76', label={'G30': 'Gaussian 530/30 nm', 'G10': 'Gaussian 530/10 nm', 'Ir(ppy)2acac': 'Ir(ppy)2acac'}[spn])
    ax[2].set_xlabel('total EQE'); ax[2].set_ylabel('on-axis current efficiency (cd/A)'); ax[2].set_yscale('log')
    ax[2].set_title('(c) on-axis cd/A vs total EQE', loc='left'); ax[2].legend(fontsize=5, frameon=False, loc='lower right')
    sh.append(('c_summary', ['spectrum', 'configuration', 'EQE_total', 'cdA_on_axis', 'du_v_30', 'du_v_60', 'half_luminance_angle'], rows))
    save(fig, 'figSel5_display', sh)
