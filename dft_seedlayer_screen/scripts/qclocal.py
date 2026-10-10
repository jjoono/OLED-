"""Local PBE0-D3(BJ)/def2 calculations with PySCF, matched to the Gaussian route.

Gaussian is the project's production code, but it runs only on the lab
workstation. This module reproduces its level -- PBE0 (PBE1PBE), def2 basis and
ECP, D3 with Becke-Johnson damping -- so single points and analyses that need no
geometry optimisation can be run anywhere. Gaussian folds the D3 energy into
"SCF Done"; here it is computed separately with s-dftd3 and added, so totals are
comparable to the Gaussian outputs line for line.
"""
import os
import numpy as np
from pyscf import gto, dft, lib
from pyscf.lo import nao as _nao, orth

BOHR = 0.52917721092
Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Mo": 42, "Ag": 47, "Li": 3,
     "Al": 13, "P": 15, "S": 16, "Cu": 29, "I": 53, "Cs": 55}


def mol(S, X, mult, basis="def2-svp", ghost=()):
    """ghost: indices whose basis functions are kept but nuclei and electrons
    removed -- the counterpoise construction."""
    atoms = []
    for i, (s, x) in enumerate(zip(S, X)):
        atoms.append((f"ghost-{s}" if i in ghost else s, tuple(x)))
    m = gto.Mole()
    m.atom = atoms
    m.unit = "Angstrom"
    m.basis = basis
    m.ecp = basis
    m.spin = mult - 1
    m.charge = 0
    m.verbose = 0
    m.max_memory = int(os.environ.get("QC_MEM_MB", "12000"))
    m.build()
    return m


def scf(m, conv=1e-9, grid=4, dm0=None):
    mf = (dft.UKS(m) if m.spin else dft.RKS(m)).density_fit()
    mf.xc = "PBE0"
    mf.grids.level = grid
    mf.conv_tol = conv
    mf.max_cycle = 200
    mf.level_shift = 0.2
    e = mf.kernel(dm0) if dm0 is not None else mf.kernel()
    if not mf.converged:
        mf = mf.newton()
        e = mf.kernel(mf.make_rdm1())
    return mf, e


def d3bj(S, X):
    from dftd3.interface import RationalDampingParam, DispersionModel
    model = DispersionModel(np.array([Z[s] for s in S]),
                            np.asarray(X, float) / BOHR)
    return model.get_dispersion(RationalDampingParam(method="pbe0"),
                                grad=False)["energy"]


def total_dm(mf):
    dm = mf.make_rdm1()
    return dm[0] + dm[1] if np.ndim(dm) == 3 else dm


def mulliken(mf):
    m = mf.mol
    dm = total_dm(mf)
    s = m.intor_symmetric("int1e_ovlp")
    pop = np.einsum("ij,ji->i", dm, s)
    q = np.array([m.atom_charge(i) for i in range(m.natm)], float)
    for i, (b0, b1) in enumerate(m.aoslice_by_atom()[:, 2:]):
        q[i] -= pop[b0:b1].sum()
    return q


def lowdin(mf, meta=False):
    m = mf.mol
    dm = total_dm(mf)
    s = m.intor_symmetric("int1e_ovlp")
    c = orth.orth_ao(m, "meta_lowdin" if meta else "lowdin", s=s)
    ci = np.linalg.solve(c, np.eye(c.shape[0]))
    p = ci @ dm @ ci.T
    pop = np.diag(p)
    q = np.array([m.atom_charge(i) for i in range(m.natm)], float)
    for i, (b0, b1) in enumerate(m.aoslice_by_atom()[:, 2:]):
        q[i] -= pop[b0:b1].sum()
    return q


def natural(mf):
    """Natural atomic orbital charges and the per-atom NAO occupations.

    Returns (charges, list of (label, occupation) per atom). This is the NPA
    of NBO: diagonalise the density in each atom's angular-momentum blocks of
    an occupancy-weighted symmetric orthogonalisation.
    """
    m = mf.mol
    dm = total_dm(mf)
    s = m.intor_symmetric("int1e_ovlp")
    c = _nao.nao(m, mf, s)
    ci = np.linalg.solve(c, np.eye(c.shape[0]))
    p = ci @ dm @ ci.T
    occ = np.diag(p)
    labels = m.ao_labels(fmt=False)
    q = np.array([m.atom_charge(i) for i in range(m.natm)], float)
    per = [[] for _ in range(m.natm)]
    for k, lab in enumerate(labels):
        q[lab[0]] -= occ[k]
        per[lab[0]].append((lab[2] + lab[3], occ[k]))
    return q, per


def iao_charges(mf):
    """Intrinsic-atomic-orbital charges (Knizia): projected onto a minimal
    free-atom basis, so nearly independent of the AO basis used."""
    from pyscf.lo import iao
    m = mf.mol
    s = m.intor_symmetric("int1e_ovlp")
    if np.ndim(mf.mo_coeff) == 3:
        occ = [mf.mo_coeff[k][:, mf.mo_occ[k] > 0] for k in (0, 1)]
        c = iao.iao(m, np.hstack(occ))
    else:
        c = iao.iao(m, mf.mo_coeff[:, mf.mo_occ > 0])
    c = orth.vec_lowdin(c, s)
    dm = total_dm(mf)
    p = c.T @ s @ dm @ s @ c
    pop = np.diag(p)
    from pyscf.lo.iao import reference_mol
    ref = reference_mol(m)
    q = np.array([m.atom_charge(i) for i in range(m.natm)], float)
    for i, (b0, b1) in enumerate(ref.aoslice_by_atom()[:, 2:]):
        q[i] -= pop[b0:b1].sum()
    return q


def ag5s(per_atom):
    """Natural occupation of the valence 5s of a silver atom: the largest
    occupied s-type NAO outside the ECP core."""
    s = sorted((o for lab, o in per_atom if lab.startswith(("4s", "5s", "s"))
                or lab[1:2] == "s"), reverse=True)
    # the ECP removes 1s-3d; the 4s,4p are semicore and stay ~2 each, so the
    # 5s is the s-type NAO whose occupation is furthest from 2
    s_like = [o for lab, o in per_atom if "s" in lab and "d" not in lab]
    cand = [o for o in s_like if o < 1.6]
    return max(cand) if cand else min(s_like)


def hirshfeld(mf, S, X, basis="def2-svp"):
    """Hirshfeld charges from spherically averaged free-atom densities at the
    same level, on the molecular grid."""
    m = mf.mol
    grids = mf.grids
    if grids.coords is None:
        grids.build()
    coords, w = grids.coords, grids.weights
    from pyscf.dft import numint
    ni = numint.NumInt()
    rho_mol = ni.get_rho(m, total_dm(mf), grids)
    free = {}
    rho_at = np.zeros((m.natm, len(w)))
    for i, s in enumerate(S):
        if s not in free:
            mult = {"H": 2, "C": 3, "N": 4, "O": 3, "F": 2, "Mo": 7,
                    "Ag": 2}.get(s, 1)
            a = mol([s], [[0, 0, 0]], mult, basis)
            af = dft.UKS(a).density_fit()
            af.xc = "PBE0"
            af.conv_tol = 1e-8
            # fractional occupation would make it spherical; a broken-symmetry
            # open shell is close enough for a Hirshfeld promolecule
            af.kernel()
            free[s] = (a, total_dm(af))
        a, dma = free[s]
        a2 = a.copy()
        a2.atom = [(s, tuple(X[i]))]
        a2.build()
        rho_at[i] = ni.get_rho(a2, dma, grids)
    pro = rho_at.sum(axis=0) + 1e-30
    q = np.array([m.atom_charge(i) for i in range(m.natm)], float)
    # ECP atoms: atom_charge already counts valence only, and the free-atom
    # densities are valence-only too, so the partition is consistent
    for i in range(m.natm):
        q[i] -= np.sum(w * rho_mol * rho_at[i] / pro)
    return q
