"""Stage 5 - deliverables: n,k tables, the summary sheet and the figures.

Everything lands in ELLIPS_OUT; the repository keeps no measurement data and no
fitted numbers.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import os, json, csv, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import jnk_data as D, jnk_ref as R, jnk_model as M, jnk_seed as S, jnk_ag as A
import jnk_tr as TR, jnk_glass as G, jnk_joint as J

OUTDIR = _os.path.join(ELLIPS_OUT, 'jnk_nk')
WLP = np.arange(300.0, 1001.0, 5.0)
HB_EVS, EPS0 = 6.582119569e-16, 8.8541878128e-12


def _load(name):
    p = _os.path.join(ELLIPS_OUT, name)
    return json.load(open(p)) if os.path.exists(p) else None


def write_nk():
    os.makedirs(OUTDIR, exist_ok=True)
    AG = _load('jnk_ag.json'); JO = _load('jnk_joint.json'); SD = _load('jnk_seed.json')
    n_files = 0
    for seed, v in SD.items():
        rec = [x for x in v['scan'] if abs(x['d_ox'] - A.D_OX) < 1e-9][0]
        N = S.seed_N(np.array(rec['p'][1:]), WLP)
        with open(_os.path.join(OUTDIR, 'seed_%s_nk.csv' % seed), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['wavelength_nm', 'n', 'k'])
            for i, wl in enumerate(WLP):
                w.writerow(['%.1f' % wl, '%.5f' % N[i].real, '%.5f' % N[i].imag])
        n_files += 1
    for sh, rec in (AG or {}).items():
        cols = {'SE': A.agN(np.array(rec['p']), WLP)}
        if JO and sh in JO:
            cols['joint'] = J.agN(np.array(JO[sh]['joint']['p']), WLP)
            cols['TR'] = J.agN(np.array(JO[sh]['tr']['p']), WLP)
        fn = 'Ag_%s_%s_%dnm_nk.csv' % (rec['seed'], sh, rec['ag_nom'])
        with open(_os.path.join(OUTDIR, fn), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['wavelength_nm'] + sum([['n_' + k, 'k_' + k] for k in cols], []))
            for i, wl in enumerate(WLP):
                row = ['%.1f' % wl]
                for k in cols:
                    row += ['%.5f' % cols[k][i].real, '%.5f' % cols[k][i].imag]
                w.writerow(row)
        n_files += 1
    print('n,k tables: %d files -> %s' % (n_files, OUTDIR))


def summary():
    AG = _load('jnk_ag.json'); JO = _load('jnk_joint.json'); TRJ = _load('jnk_tr.json')
    if not AG:
        return
    ref = max((v['wp'] for v in AG.values()), default=1.0)
    print('\n%-6s %-6s %4s | %6s %6s %6s | %6s %6s | %6s %6s | %5s %5s'
          % ('sheet', 'seed', 'nom', 'd_SE', 'd_TR', 'd_join', 'n633', 'k633',
             'hw_p', 'rho_f', 'MSE', 'dT%p'))
    for sh, rec in AG.items():
        j = (JO or {}).get(sh, {})
        t = (TRJ or {}).get(sh, {})
        dj = j.get('joint', {}).get('p', [np.nan])[0]
        print('%-6s %-6s %4d | %6.2f %6.2f %6.2f | %6.3f %6.3f | %6.2f %6.2f | %5.2f %5.2f'
              % (sh, rec['seed'], rec['ag_nom'], rec['p'][0], t.get('d_TR', np.nan), dj,
                 rec['n'][2], rec['k'][2], rec['wp'], (rec['wp'] / ref)**2,
                 rec['mse'], j.get('joint', {}).get('dT', np.nan)))


NKLIB = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                      '..', '..', '..', 'dft_seedlayer_screen', 'data', 'nk')


def cross_check():
    """Against the optical constants the repository already carries."""
    AG = _load('jnk_ag.json'); JO = _load('jnk_joint.json'); SD = _load('jnk_seed.json')
    if not AG or not os.path.isdir(NKLIB):
        return
    lib = lambda f: np.genfromtxt(_os.path.join(NKLIB, f), delimiter=',', names=True)
    print('\nAg on HATCN: this work vs data/nk/Ag<d>nm_on_HATCN5_measured.csv '
          '(the campaign T/R inversion)')
    print('%4s %6s | %-24s | %-24s'
          % ('Ag', 'wl', 'n   SE / joint / libTR', 'k   SE / joint / libTR'))
    for nom, sh in ((5, '1-7'), (7, '2-5'), (8, '2-6')):
        fn = 'Ag%dnm_on_HATCN5_measured.csv' % nom
        if not os.path.exists(_os.path.join(NKLIB, fn)) or sh not in AG:
            continue
        d = lib(fn)
        for w in (450., 550., 650., 750.):
            Nse = A.agN(np.array(AG[sh]['p']), np.array([w]))[0]
            Njo = (J.agN(np.array(JO[sh]['joint']['p']), np.array([w]))[0]
                   if JO and sh in JO else complex(np.nan, np.nan))
            nl = np.interp(w, d['wavelength_nm'], d['n'])
            kl = np.interp(w, d['wavelength_nm'], d['k'])
            print('%4d %6.0f | %7.3f %7.3f %7.3f  | %7.3f %7.3f %7.3f'
                  % (nom, w, Nse.real, Njo.real, nl, Nse.imag, Njo.imag, kl))
    if SD and os.path.exists(_os.path.join(NKLIB, 'l_HATCN.csv')):
        d = lib('l_HATCN.csv')
        rec = [x for x in SD['HATCN']['scan'] if abs(x['d_ox'] - A.D_OX) < 1e-9][0]
        print('\nHATCN seed: this SE fit vs data/nk/l_HATCN.csv')
        for w in (450., 550., 633., 750.):
            N = S.seed_N(np.array(rec['p'][1:]), np.array([w]))[0]
            print('  %5.0f nm   n %7.4f / %7.4f     k %8.5f / %8.5f'
                  % (w, N.real, np.interp(w, d['wavelength_nm'], d['n']),
                     N.imag, np.interp(w, d['wavelength_nm'], d['k'])))


def figures():
    AG = _load('jnk_ag.json'); JO = _load('jnk_joint.json'); SD = _load('jnk_seed.json')
    if not AG:
        return
    os.makedirs(OUTDIR, exist_ok=True)
    Nb = R.n_ag_bulk(WLP)

    fig, ax = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    for col, seed in enumerate(('HATCN', 'MoOx')):
        rows = [(sh, v) for sh, v in AG.items() if v['seed'] == seed]
        rows.sort(key=lambda x: x[1]['ag_nom'])
        cm = plt.cm.viridis(np.linspace(0, 0.9, len(rows)))
        for c, (sh, v) in zip(cm, rows):
            if JO and sh in JO:
                N = J.agN(np.array(JO[sh]['joint']['p']), WLP)
            else:
                N = A.agN(np.array(v['p']), WLP)
            ax[0, col].plot(WLP, N.real, color=c, label='%d nm' % v['ag_nom'])
            ax[1, col].plot(WLP, N.imag, color=c)
        for r in (0, 1):
            ax[r, col].plot(WLP, Nb.real if r == 0 else Nb.imag, 'k--', lw=1, label='bulk Ag')
            ax[r, col].grid(alpha=.3)
        ax[0, col].set_title('Ag on %s  (joint SE + T/R fit)' % seed)
        ax[0, col].set_ylabel('n'); ax[1, col].set_ylabel('k')
        ax[1, col].set_xlabel('wavelength (nm)')
        ax[0, col].legend(fontsize=7, ncol=2)
    fig.tight_layout(); fig.savefig(_os.path.join(OUTDIR, 'nk_by_seed.png'), dpi=140)
    plt.close(fig)

    if JO:
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for col, seed in enumerate(('HATCN', 'MoOx')):
            rows = sorted([v for v in JO.values() if v['seed'] == seed], key=lambda x: x['ag_nom'])
            nom = [v['ag_nom'] for v in rows]
            for key, mk in (('se', 'o-'), ('joint', 's-'), ('tr', '^-')):
                ax[col].plot(nom, [v[key]['n'][1] for v in rows], mk, label='n(550) %s' % key)
            ax[col].set_title('Ag on %s' % seed)
            ax[col].set_xlabel('nominal Ag (nm)'); ax[col].set_ylabel('n at 550 nm')
            ax[col].grid(alpha=.3); ax[col].legend(fontsize=8)
        fig.tight_layout(); fig.savefig(_os.path.join(OUTDIR, 'n550_se_vs_tr.png'), dpi=140)
        plt.close(fig)

    if SD:
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for seed, v in SD.items():
            rec = [x for x in v['scan'] if abs(x['d_ox'] - A.D_OX) < 1e-9][0]
            N = S.seed_N(np.array(rec['p'][1:]), WLP)
            ax[0].plot(WLP, N.real, label='%s (%.1f nm)' % (seed, rec['p'][0]))
            ax[1].plot(WLP, N.imag, label=seed)
        for a, lab in zip(ax, ('n', 'k')):
            a.set_xlabel('wavelength (nm)'); a.set_ylabel(lab); a.grid(alpha=.3); a.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(_os.path.join(OUTDIR, 'seed_nk.png'), dpi=140)
        plt.close(fig)
    if JO:
        import jnk_tr as TRm
        fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
        for col, seed in enumerate(('HATCN', 'MoOx')):
            rows = sorted([v for v in JO.values() if v['seed'] == seed], key=lambda x: x['ag_nom'])
            cm = plt.cm.viridis(np.linspace(0, 0.9, len(rows)))
            for c, v in zip(cm, rows):
                wl, Tm, Rm = TRm.measured(v['tag'])
                p = np.array(v['joint']['p'])
                Na = J.agN(p, wl)
                Ns = S.seed_N(A.seed_params()[seed][1], wl)
                layers, ds = TRm.stack_of(Na, p[2], p[3], Ns, v['d_seed'], wl)
                T, Rr = M.stack_TR(wl, layers, ds, G.n_glass(wl), tau=G.tau_glass(wl))
                if Tm is not None:
                    ax[0, col].plot(wl, 100 * Tm, color=c, lw=2.4, alpha=.35)
                    ax[0, col].plot(wl, 100 * T, color=c, lw=1, label='%d nm' % v['ag_nom'])
                if Rm is not None:
                    ax[1, col].plot(wl, 100 * Rm, color=c, lw=2.4, alpha=.35)
                    ax[1, col].plot(wl, 100 * Rr, color=c, lw=1)
            ax[0, col].set_title('Ag on %s  (thick = measured, thin = joint model)' % seed, fontsize=10)
            ax[0, col].set_ylabel('T (%)'); ax[1, col].set_ylabel('R (%)')
            ax[1, col].set_xlabel('wavelength (nm)')
            for r in (0, 1):
                ax[r, col].grid(alpha=.3)
            ax[0, col].legend(fontsize=7, ncol=2)
        fig.tight_layout(); fig.savefig(_os.path.join(OUTDIR, 'tr_joint_fit.png'), dpi=140)
        plt.close(fig)
    print('figures -> %s' % OUTDIR)


if __name__ == '__main__':
    write_nk()
    summary()
    cross_check()
    figures()
