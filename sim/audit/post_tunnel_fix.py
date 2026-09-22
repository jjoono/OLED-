"""After the substrate-tunnelling fix in cps2.solve_pol: print what moved, rebuild the workbooks and Table 1."""
import numpy as np, subprocess, re, csv
R = '/home/user/OLED-'
D = np.genfromtxt(R + '/sim/audit/fig2def_series.csv', delimiter=',', names=True, dtype=None, encoding='utf-8')
print('Fig 2(d)-(f) rows with n_sub >= 1.80 (eta_sub, spp, abs, Aprime, eta_ext_closed, EQE_closed):')
for r in D:
    if r['n_sub'] >= 1.795: print('  %s %.2f  %.4f %.4f %.4f %.4f  %.4f %.4f' % (r['reflector'], r['n_sub'], r['eta_sub'], r['spp'], r['abs'], r['Aprime'], r['eta_ext_closed'], r['EQE_closed']))
T = list(csv.DictReader(open(R + '/sim/audit/table1_series.csv')))
print('Table 1 (eta_sub, spp, Aprime, eta_ext_series, EQE_series, eta_ext_closed, EQE_closed):')
for r in T: print('  %s %s %s  %.3f %.3f %.3f  %.3f %.3f  %.3f %.3f' % (r['d_ETL_nm'], r['k_ITO'], r['Theta'], float(r['eta_sub']), float(r['spp']), float(r['Aprime']), float(r['eta_ext_series']), float(r['EQE_series']), float(r['eta_ext_closed']), float(r['EQE_closed'])))
C = np.genfromtxt(R + '/sim/audit/fig5c_series.csv', delimiter=',', names=True, dtype=None, encoding='utf-8')
print('Fig 5(c) proposed at Theta 0.5/0.65/0.70/1.0:', ['%.4f' % r['EQE_series'] for r in C if r['device'] == 'proposed' and r['Theta'] in (0.5, 0.65, 0.7, 1.0)])
subprocess.run(['python3', 'make_update_xlsx.py'], cwd=R + '/figures/update', check=True)
subprocess.run(['python3', 'make_audit_xlsx.py'], cwd=R + '/sim/audit', check=True)
subprocess.run(['python3', 'make_table1.py'], cwd=R + '/figures/table1', check=True)
# Table 1 rows in the manuscript
p = R + '/manuscript/unityEQE_final_ko.md'; s = open(p, encoding='utf-8').read()
i0 = s.index('| d_ETL (nm) | k_ITO |'); i1 = s.index('\n\n', i0)
rows = ['| d_ETL (nm) | k_ITO | Θ | η_sub | SPP 손실 | A′ | η_ext | EQE |', '|---|---|---|---|---|---|---|---|']
for r in T:
    th = float(r['Theta']); rows.append('| %s | %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f |' % (r['d_ETL_nm'], r['k_ITO'], '0.67' if abs(th - 0.667) < 1e-3 else '%.2f' % th, float(r['eta_sub']), float(r['spp']), float(r['Aprime']), float(r['eta_ext_series']), float(r['EQE_series'])))
s = s[:i0] + '\n'.join(rows) + s[i1:]
ec = [float(r['eta_ext_closed']) for r in T]; qc = [float(r['EQE_closed']) for r in T]; es = [float(r['eta_ext_series']) for r in T]; qs = [float(r['EQE_series']) for r in T]
foot = '식 (3)으로 계산하면 η_ext %.3f–%.3f, EQE %.3f–%.3f로 급수보다 %.1f–%.1f%%p 높다(식 (3)은 상한).' % (min(ec), max(ec), min(qc), max(qc), 100 * min(np.array(qc) - np.array(qs)), 100 * max(np.array(qc) - np.array(qs)))
s = re.sub(r'식 \(3\)으로 계산하면 η_ext [^.]*\.', foot, s, count=1)
open(p, 'w', encoding='utf-8').write(s); print('manuscript table refreshed; footnote:', foot)
