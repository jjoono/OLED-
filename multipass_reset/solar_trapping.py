"""Light trapping check: does the placement of angular mixing matter in a weakly absorbing slab (solar-cell geometry)?

Slab of index n with a textured front face (traced matrices from run_highindex.py) and a back reflector of
reflectance Rb.  Light arrives from air at normal incidence.  Absorption per one-way crossing at internal angle
theta is 1 - exp(-alpha d / cos theta).  Compared (same front-face transmission p(theta) in every case):
  A  actual texture, specular back reflector
  B  front face mixes its returned light completely in reflection (rank-one, q ~ w(1-p)), specular back
  C  actual texture, Lambertian back reflector (mixing at the back / in transmission)
  D  front face reflective mixing + Lambertian back
  L  Lambertian reference (Yablonovitch): weak-absorption path enhancement 4 n^2
Path enhancement F = absorbed fraction / (alpha d) in the weak-absorption limit.
Output: ../selectivity/data/solar_trapping.json (gitignored) and a printed table."""
import json, numpy as np

E = np.linspace(0, 90, 91); TH = np.radians(0.5 * (E[:-1] + E[1:]))
W = np.diff(np.sin(np.radians(E)) ** 2)                          # etendue per bin (per unit n^2)
COS = np.cos(TH)
CONE = 10


def load(name, n):
    z = np.load(f'highindex/{name}_n{n}.npz')
    return z['BT'], z['BR'].T                                   # BT[in, air]; BR[out, in]


def planar(n):
    """Planar face: Fresnel (unpolarized), specular.  Internal bin k maps to the air bins it overlaps in
    n^2 sin^2(theta) (etendue), so the binned matrix keeps reciprocity w_air T(air->in) = n^2 w_in T(in->air)."""
    s = n * np.sin(TH); ok = s < 1; ct = np.sqrt(np.clip(1 - s ** 2, 0, None)); ci = COS
    rs = ((n * ci - ct) / (n * ci + ct)) ** 2; rp = ((ci - n * ct) / (ci + n * ct)) ** 2
    T = np.where(ok, 1 - 0.5 * (rs + rp), 0.0)
    S2 = np.sin(np.radians(E)) ** 2; BT = np.zeros((90, 90))
    for k in np.where(ok)[0]:
        lo, hi = n * n * S2[k], min(n * n * S2[k + 1], 1.0)
        ov = np.clip(np.minimum(S2[1:], hi) - np.maximum(S2[:-1], lo), 0, None)
        if ov.sum() > 0: BT[k] = T[k] * ov / ov.sum()
    return BT, np.diag(1 - T)


def run(BT, BR, n, ad, Rb=1.0, front='actual', back='specular'):
    p = BT.sum(1)
    # near-normal sunlight (air cone 0-10 deg, uniform in etendue) coupled into internal bin k by reciprocity:
    # w_air,j T(air j -> in k) = n^2 w_k BT[k, j]; normalized per unit incident power
    din = n ** 2 * W * BT[:, :CONE].sum(1) / W[:CONE].sum()
    if front == 'mix':
        q = W * (1 - p); q = q / q.sum(); BR = (1 - p)[None, :] * q[:, None]
    one = np.exp(-ad / COS)
    absorbed = 0.0; d = din.copy(); escaped = 0.0
    for _ in range(20000):
        # down through the slab, back reflector, up through the slab
        absorbed += (d * (1 - one)).sum(); d = d * one
        if back == 'lambert':
            tot = Rb * d.sum(); u = tot * W / W.sum()
        else:
            u = Rb * d
        absorbed += (u * (1 - one)).sum(); u = u * one
        escaped += (u * p).sum(); d = BR @ u
        if d.sum() < 1e-12: break
    return absorbed, escaped, din.sum()


out = {}
for n in (1.51, 2.4, 3.5):
    for name in ('planar', 'hemiMLA', 'wavy0.4', 'wavy1.0'):
        BT, BR = planar(n) if name == 'planar' else load(name, n)
        p = BT.sum(1); wn = W / W.sum(); pm = p @ wn; s2 = ((p - pm) ** 2) @ wn; pmix = pm - s2 / (1 - pm)
        row = dict(p=float(pm), p_mix=float(pmix), sigma2=float(s2), coupled=float(n ** 2 * (W * BT[:, :CONE].sum(1)).sum() / W[:CONE].sum()))
        for ad in (1e-3, 0.05):
            for key, fr, bk in (('A', 'actual', 'specular'), ('B', 'mix', 'specular'), ('C', 'actual', 'lambert'), ('D', 'mix', 'lambert')):
                a, e, c = run(BT, BR, n, ad, front=fr, back=bk)
                row[f'{key}_ad{ad:g}'] = dict(absorbed=a, F=a / ad, escaped=e)
        out[f'{name}|n{n}'] = row
        r = row
        print(f"n={n:4} {name:8s} p {r['p']:.3f} p_mix {r['p_mix']:.3f} coupled {r['coupled']:.3f} | F (ad=1e-3): "
              f"A {r['A_ad0.001']['F']:6.1f}  B front-mix {r['B_ad0.001']['F']:6.1f}  C back-Lamb {r['C_ad0.001']['F']:6.1f}  D both {r['D_ad0.001']['F']:6.1f}"
              f"  | 4n^2 = {4 * n * n:.1f} | absorbed (ad=0.05): A {r['A_ad0.05']['absorbed']:.3f} B {r['B_ad0.05']['absorbed']:.3f} C {r['C_ad0.05']['absorbed']:.3f} D {r['D_ad0.05']['absorbed']:.3f}", flush=True)
json.dump(out, open('../selectivity/data/solar_trapping.json', 'w'), indent=1)
