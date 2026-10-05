"""Fig. 1b,c raw data with round-trip losses tied to Fig. 2 (A' = 0.02 ideal, 0.045 Ag, 0.15 Al; ITO 150 nm, Fig. 2g-j)."""
import os, numpy as np, openpyxl
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
P = 0.4
CASES = [(0.02, 'low-loss limit'), (0.045, 'Ag reflector (Fig. 2j)'), (0.15, 'Al reflector (Fig. 2j)')]
K = np.arange(0, 21)
eta = lambda p, a: p / (p + (1 - p) * a)
wb = openpyxl.Workbook(); ws = wb.active; ws.title = 'README'
for l in ["Fig. 1b,c raw data. Eq. (2): eta_sta = p / (p + (1 - p) A').",
          "Fraction extracted in round trip k (k = 0 is the first encounter): p [(1 - p)(1 - A')]^k; cumulative = sum over passes up to k.",
          "A' = 0.02: low-loss limit; 0.045: Ag reflector with 150 nm ITO; 0.15: Al reflector with 150 nm ITO (spectrum- and cos(theta)sin(theta)-weighted round-trip losses 4.5% and 15.6% of Fig. 2j).",
          "p = 0.4 (glass, ~1/n^2). Panel c: eta_sta vs R_LED = 1 - A' for p = 0.25, 0.4, 0.6."]:
    ws.append([l])
b = wb.create_sheet('b_per_round_trip')
b.append(['round trip k'] + [f"A'={a} per pass" for a, _ in CASES] + [f"A'={a} cumulative" for a, _ in CASES])
for k in K:
    per = [P * ((1 - P) * (1 - a)) ** k for a, _ in CASES]
    cum = [sum(P * ((1 - P) * (1 - a)) ** j for j in range(k + 1)) for a, _ in CASES]
    b.append([int(k)] + [round(x, 6) for x in per] + [round(x, 6) for x in cum])
b.append([]); b.append(['limit (Eq. 2)', '', '', ''] + [round(eta(P, a), 4) for a, _ in CASES])
c = wb.create_sheet('c_eta_vs_R')
R = np.round(np.linspace(0.70, 1.0, 301), 4)
c.append(['R_LED', 'p = 0.25', 'p = 0.4', 'p = 0.6'])
for r in R: c.append([float(r)] + [round(eta(p, 1 - r), 6) for p in (0.25, 0.4, 0.6)])
c.append([]); c.append(['marked points (p = 0.4)', 'R_LED', 'eta_sta'])
for a, lab in CASES: c.append([lab, 1 - a, round(eta(P, a), 4)])
wb.save(os.path.join(HERE, 'fig1bc_rawdata_v2.xlsx'))
fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
cols = ['#0072B2', '#56B4E9', '#C8553D']
for (a, lab), col, off in zip(CASES, cols, (-0.27, 0, 0.27)):
    per = P * ((1 - P) * (1 - a)) ** K
    ax[0].bar(K[:11] + off, per[:11], 0.27, color=col, alpha=0.6)
    ax[0].plot(K[:11], np.cumsum(per)[:11], color=col, label=f"A' = {a} ({eta(P,a):.2f})")
for p, ls in ((0.25, '--'), (0.4, '-'), (0.6, ':')):
    ax[1].plot(R, eta(p, 1 - R), 'k', ls=ls, lw=1, label=f'p = {p}')
for (a, lab), col in zip(CASES, cols): ax[1].plot(1 - a, eta(P, a), 'o', color=col)
ax[0].set(xlabel='round trip', ylabel='extracted fraction', ylim=(0, 1)); ax[0].legend(fontsize=7, frameon=False)
ax[1].set(xlabel='$R_{LED}$ = 1 − A′', ylabel=r'$\eta_{sta}$', xlim=(0.7, 1), ylim=(0.2, 1)); ax[1].legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(os.path.join(HERE, 'fig1bc_v2_preview.png'), dpi=200)
print({a: round(eta(P, a), 3) for a, _ in CASES})
