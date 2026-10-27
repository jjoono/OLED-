"""Expected s-polarised reflectance at 30 deg (incidence from air on the glass side, n 1.52) of the reference-device pixel,
for the absorption scenarios of absorb_scan.py (stack refitted to the angular EL in each case).
Stack-only R (inside glass, 19.2 deg) and the measured total incl. the air/glass front face (incoherent).  -> refl_30s.csv"""
import sys, csv, numpy as np
sys.path[:0] = ['/home/user/OLED-/sim/diag', '/home/user/OLED-/sim/fig3c']
import fit_ref_angles as F, cps2
NS, A = 1.52, 30.0
TG = np.degrees(np.arcsin(np.sin(np.radians(A)) / NS))
c, ct = np.cos(np.radians(A)), np.cos(np.radians(TG)); RF = ((c - NS * ct) / (c + NS * ct)) ** 2      # s, air/glass
CASES = {  # name: (f, g, z, IZO, kd, kh, kcs)   from absorb_scan(_cs).csv fits; nominal = as designed
    'nominal (ETL185, IZO50)': (1.0, 1.0, 1.0, 50.0, 0, 0, 0),
    'no extra absorption':      (172.5 / 185, 160.0 / 210, 1.0, 12.7, 0, 0, 0),
    'kd 0.02 (doped ETL)':      (172.5 / 185, 156.9 / 210, 1.0, 14.5, 0.02, 0, 0),
    'kd 0.04 (doped ETL)':      (172.4 / 185, 154.1 / 210, 1.0, 16.3, 0.04, 0, 0),
    'kh 0.1 (HATCN/TAPC)':      (172.5 / 185, 156.4 / 210, 1.0, 14.5, 0, 0.1, 0),
    'kh 0.2 (HATCN/TAPC)':      (172.5 / 185, 153.5 / 210, 1.0, 16.1, 0, 0.2, 0),
    'kcs 0.6 (Cs2CO3/Ag)':      (172.5 / 185, 156.5 / 210, 1.0, 14.6, 0, 0, 0.6),
}
I = np.clip(np.asarray(F.D.spec.l_I_Irppy2acac, float)[:301], 0, None); W = I * F.D.LAM; W /= W.sum()
rows = []
for name, (f, g, z, di, kd, kh, kc) in CASES.items():
    Rs = []
    for j in range(301):
        S = F.stack(j, f, g, z, di, 1.0, kd, kh, kc); k0 = 2 * np.pi / S.lam
        kx = NS * np.sin(np.radians(TG)) * k0
        rs, _ = cps2.stack_rt(S.all_no, S.all_ne, S.all_d, np.array([kx]), k0, 's'); Rs.append(float(np.abs(rs[0]) ** 2))
    Rs = np.array(Rs); Rt = RF + (1 - RF) ** 2 * Rs / (1 - RF * Rs)
    r = dict(case=name, stack_avg=round(float(W @ Rs), 4), total_avg=round(float(W @ Rt), 4),
             **{'total_%d' % l: round(float(Rt[l - 400]), 4) for l in (450, 500, 520, 550, 600, 650, 700)},
             total_min=round(float(Rt.min()), 4), lam_min=int(400 + Rt.argmin()))
    rows.append(r); print(r, flush=True)
    np.savetxt('refl_30s_%s.txt' % name.split()[0].replace('(', ''), np.c_[F.D.LAM, Rs, Rt], fmt='%.5f', header='lam R_stack_s R_total_s')
with open('refl_30s.csv', 'w', newline='') as fo:
    w = csv.DictWriter(fo, list(rows[0])); w.writeheader(); w.writerows(rows)
print('front-face s reflectance %.4f, angle in glass %.1f deg' % (RF, TG))
