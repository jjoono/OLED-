"""EQE of the user's BE device on an n 1.77 substrate with the cone array (infinite, thin, lossless substrate; matrix series).
  air | Ag 100 | B3PyMPM 170 | B3PyMPM 10 | TCTA:B3PyMPM:Ir(ppy)2acac 25 (dipole at centre) | TCTA 10 | TAPC 200 | HATCN 5
  | ITO 50 (l_ITO_SNU_temp) | substrate n 1.77 | cone array (r 1.4 um, h 1.8 um, 0.2 um gap) | air
Ir(ppy)2acac 400-700 nm, photon-number weighted, 2 nm grid; PLQY 1 (scale linearly)."""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import cps2, series
M = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)
mat, spec = M['material'], M['spectrum']
LAM = np.arange(400, 701)
def lib(f): return np.asarray(getattr(mat, f))[:301]
AG, ITO, HAT = lib('l_Ag_McPeak'), lib('l_ITO_SNU_temp'), lib('l_HATCN')
B3o, B3e = lib('l_B3_o_JO'), lib('l_B3_e_JO'); EMo, EMe = lib('l_TCTA_B3_o_JO'), lib('l_TCTA_B3_e_JO')
TCo, TCe = lib('l_TCTA_o_JO'), lib('l_TCTA_e_JO'); TPo, TPe = lib('l_TAPC_o_JO'), lib('l_TAPC_e_JO')
I = np.clip(np.asarray(spec.I_Irppy2acac, float)[:301], 0, None); WP = LAM * I
NS = 1.77; TH = np.arange(90) + 0.5
def split(B): return B[:90].sum(0), B[90:180][::-1]
C = np.load('/home/user/OLED-/sim/mla/cone_bsdf.npz')
c = np.cos(np.radians(TH)); s2 = (NS * np.sin(np.radians(TH))) ** 2; ok = s2 < 1; ct = np.sqrt(np.clip(1 - s2, 0, 1))
Rf = np.where(ok, 0.5 * (((NS * c - ct) / (NS * c + ct)) ** 2 + ((NS * ct - c) / (NS * ct + c)) ** 2), 1.0)
EXT = {'cone': split(C['cone_n1.77']), 'hemisphere (close-packed)': split(C['hemi_cp_n1.77']), 'flat (no extractor)': (1 - Rf, np.diag(Rf))}


def one(a):
    h, j = a
    S = cps2.Stack(float(LAM[j]), (EMo[j], EMe[j]), 25.0, 12.5,
                   above=[(B3o[j], B3e[j], 10.0), (B3o[j], B3e[j], 170.0), (AG[j], AG[j], 100.0)],
                   below=[(TCo[j], TCe[j], 10.0), (TPo[j], TPe[j], 200.0), (HAT[j], HAT[j], 5.0), (ITO[j], ITO[j], 50.0)],
                   n_sub=NS, h=h)
    r = cps2.solve_pol(S, npts=3000); R = np.real(cps2.stack_reflectance(S, TH)); P = np.real(cps2.sub_angular(S, TH))
    es = float(np.real(r['air'] + r['sub'])); Pn = P / P.sum(); Ap = float(1 - Pn @ R)
    return h, j, es, float(np.real(r['spp'])), float(np.real(r['wg'])), float(np.real(r['abs'])), Ap, [series.eta_ext(BT, BR, R, P, 300)[0] for BT, BR in EXT.values()]


if __name__ == '__main__':
    js = list(range(0, 301, 2)); w = WP[js] / WP[js].sum()
    thetas = (2 / 3, 0.77, 0.8)
    with Pool(4) as p:
        res = p.map(one, [(h, j) for h in thetas for j in js])
    rows = []
    for h in thetas:
        blk = [r for r in res if r[0] == h]
        es = np.array([b[2] for b in blk]); ee = np.array([b[7] for b in blk])
        sp, wg, ab, Ap = (w @ np.array([b[k] for b in blk]) for k in (3, 4, 5, 6))
        r550 = [b for b in blk if LAM[b[1]] == 550][0]
        for i, name in enumerate(EXT):
            eqe = float(w @ (es * ee[:, i]))
            rows.append(dict(Theta=round(h, 3), extractor=name, eta_sub=float(w @ es), spp=sp, wg=wg, abs_dev=ab, Aprime=Ap,
                             eta_ext=eqe / float(w @ es), EQE=eqe, EQE_550=r550[2] * r550[7][i]))
            print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
    with open('device_eqe.csv', 'w', newline='') as f:
        wr = csv.DictWriter(f, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
