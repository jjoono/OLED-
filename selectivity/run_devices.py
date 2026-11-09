"""Device optics (eta_sub, P_sub, R_LED per wavelength) for every stack, and the slab randomizers."""
import numpy as np, time, sel, slab

if __name__ == '__main__':
    t0 = time.time(); out = {}
    for k in sel.STACKS:
        d = sel.device(k)
        for f in ('eta_sub', 'P', 'R'):
            out[f'{k}_{f}'] = d[f]
        out[f'{k}_n_sub'] = d['n_sub']
        A = 1 - (d['R'] * sel.WL).sum(1)
        print(f"{k:8s} n_sub {d['n_sub']}  eta_sub(530) {d['eta_sub'][65]:.3f}  A'(Lambertian, 530 nm) {A[65]:.4f}  [{time.time() - t0:.0f}s]", flush=True)
    np.savez('data/devices.npz', **out)
    sl = {}
    for g in (0.0, 0.5, 0.8, 0.9):
        for tau in (0.1, 0.25, 0.5, 1.0, 2.0, 4.0):
            S = tau / (1 - g)
            Ts, Rs = slab.slab(S, g, 1.0, N=20000)
            sl[f'g{g}_t{tau}_w1'] = np.stack([Ts, Rs])
    for w in (0.9995, 0.999, 0.995):          # absorbing host/particles, at the reference slab
        Ts, Rs = slab.slab(5.0, 0.8, w, N=20000)
        sl[f'g0.8_t1.0_w{w}'] = np.stack([Ts, Rs])
    np.savez('data/slabs.npz', **sl)
    print('slabs done', len(sl), f'[{time.time() - t0:.0f}s]')
