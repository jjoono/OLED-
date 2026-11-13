"""Real scattering layers at the escape surface versus at the device side (ray optics, adding method).

The six non-absorbing Henyey-Greenstein slabs of the main text are placed either
  exit side  : between the substrate and the planar glass/air face (as in the main-text set), or
  device side: between the OLED and the substrate, with a planar glass/air face (and, separately, under the
               hemispherical and a random MLA).
Device side: emission enters the slab from the OLED, v_eff = Ts (I - diag R Rs)^-1 P; light returning from the
substrate meets D_eff = Rs + Ts diag R (I - Rs diag R)^-1 Ts, which replaces diag R in the multi-pass series.
All calculations at 0.5 deg (planar face) or 1 deg (MLA) substrate bins; Al, Ag, low-loss stacks at 550 nm.
Writes data/device_side_slabs.json."""
import json, numpy as np
import sel, mixing
src = open('oe_figs.py').read()
exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])

NB = sel.NB; I2 = np.eye(NB); I1 = np.eye(90)
Rf, Tf = sel.tmm_RT([], 550.0, sel.TH, N_SUB)
w2 = sel.WL / sel.WL.sum(); W2 = sel.WL.reshape(90, 2); fpair = (W2 / W2.sum(1, keepdims=True)).ravel()
def coarse_op(B):   # [out, in] 0.5 deg -> 1 deg, Lambertian within each input bin
    return (B * fpair[None, :]).reshape(90, 2, 90, 2).sum(axis=(1, 3))

out = {}
for dev in DEVS:
    i = [k for k, lam in WLS if lam == 550.0][0]
    P = DV[f'{dev}_P'][i]; R = DV[f'{dev}_R'][i]            # 0.5 deg bins
    # no scatterer, planar face
    x = np.linalg.solve(I2 - np.diag(R * Rf), P); none = float(Tf @ x)
    res = {'none': none}
    for k in ('g0.0_t0.25_w1', 'g0.0_t1.0_w1', 'g0.8_t0.25_w1', 'g0.8_t1.0_w1', 'g0.9_t1.0_w1', 'g0.8_t4.0_w1'):
        Ts, Rs = SL[k]
        # exit side: slab under the planar face
        U = np.linalg.solve(I2 - Rs * Rf[None, :], Ts); BR = Rs + (Ts * Rf[None, :]) @ U
        xe = np.linalg.solve(I2 - R[:, None] * BR, P); exit_ = float(Tf @ (U @ xe))
        # device side: slab on the OLED, planar face
        Ud = np.linalg.solve(I2 - Rs * R[None, :], Ts); Deff = Rs + (Ts * R[None, :]) @ Ud
        veff = Ts @ np.linalg.solve(I2 - R[:, None] * Rs, P)          # up through the slab after bounces between slab and OLED
        xd = np.linalg.solve(I2 - Deff * Rf[None, :], veff); devside = float(Tf @ xd)
        # device side under MLAs (1 deg)
        Dc = coarse_op(Deff); vc = veff.reshape(90, 2).sum(1)
        mla = {}
        for nm in ('MLA hemi_total', 'A_shape0.7'):
            c, BT, BRm = LIB[nm]; pT = BT.sum(1)
            Rc = R.reshape(90, 2).mean(1); Pc = P.reshape(90, 2).sum(1)
            alone = float(pT @ np.linalg.solve(I1 - Rc[:, None] * BRm, Pc))
            withs = float(pT @ np.linalg.solve(I1 - Dc @ BRm, vc))
            mla[nm] = dict(alone=alone, with_device_slab=withs)
        res[k.replace('_w1', '')] = dict(exit_side=exit_, device_side=devside, mla=mla)
    out[dev] = res
json.dump(out, open('data/device_side_slabs.json', 'w'), indent=1)
for dev, r in out.items():
    print(f"{dev}: planar, no scatterer {r['none']:.3f}")
    for k, v in r.items():
        if k == 'none': continue
        m = v['mla']
        print(f"  {k:10s} exit side {v['exit_side']:.3f} | device side {v['device_side']:.3f} | hemi MLA {m['MLA hemi_total']['alone']:.3f} -> {m['MLA hemi_total']['with_device_slab']:.3f} | random MLA {m['A_shape0.7']['alone']:.3f} -> {m['A_shape0.7']['with_device_slab']:.3f}")
