"""Kernels of the green-device freeform lenslet (cup with central spike and undercut; profile from the STL, radius
10.2 um, height 14 um, flat glass between lenslets) and validation of revtrace on the hemisphere against hftrace."""
import sys, time, numpy as np
sys.path.insert(0, '../multipass_reset'); import revtrace as RT
w = np.diff(np.sin(np.radians(np.arange(91))) ** 2)
what = sys.argv[1] if len(sys.argv) > 1 else 'both'
if what in ('both', 'validate'):
    a = np.linspace(0, np.pi / 2, 400); t = time.time()
    BT, K, lost, KT = RT.kernel(np.sin(a), np.cos(a), N=20000)
    ref = np.load('kernel_hemisphere.npz')['BT'].sum(1)
    print(f'hemisphere via revtrace: {time.time()-t:.0f}s p {np.sum(BT.sum(1)*w):.4f} (hftrace {np.sum(ref*w):.4f}), max |dp| {np.abs(BT.sum(1)-ref).max():.3f}, lost {lost.max():.4f}', flush=True)
if what in ('both', 'green'):
    P = np.load(sys.argv[2] if len(sys.argv) > 2 else 'green_profile_rz.npy')
    P = np.r_[P, [[10.2, 0.0]]] / 10.2                     # close the solid with the outer cylindrical wall
    t = time.time(); BT, K, lost, KT = RT.kernel(P[:, 0], P[:, 1], N=20000)
    np.savez('kernel_green_freeform.npz', BT=BT, K=K, lost=lost, KT=KT)
    print(f'green freeform: {time.time()-t:.0f}s p {np.sum(BT.sum(1)*w):.4f} lost max {lost.max():.4f}', flush=True)
