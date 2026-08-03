"""Build periodic slab models and write CP2K / Quantum ESPRESSO inputs.

WHY: every DFT number in this project so far comes from a molecular CLUSTER
model (REPORT.md sec.1, "한계"). Clusters cannot give a diffusion barrier on an
extended surface -- the TPBi and B3PyMPM barriers failed outright because the
bridge site has no meaning in a single molecule -- and they leave the absolute
adsorption energies on ionic solids uncertain. Periodic slabs fix both.

SYSTEMS
  hatcn_ml      HATCN monolayer, flat-lying, hexagonal 2D cell + vacuum
  hatcn_ml_ag   ... + one Ag adatom on a nitrile N (the E_b = 1.03 eV site)
  hatcn_ml_ag_ring  ... + Ag over the aromatic core (the ~0 eV site)
  hatcn_neb     7-image path N-site -> bridge -> N-site, for CI-NEB
  lif001        LiF(001) 2x2 slab, 4 layers  (cheap workflow validation)
  lif001_ag     ... + Ag adatom on F-top

CODE CHOICE IN THIS ENVIRONMENT
  CP2K runs here: the GTH pseudopotentials and MOLOPT basis sets it needs ship
  inside the pyscf wheel (pyscf/pbc/gto/{pseudo,basis}) in exactly CP2K's format,
  and PyPI is reachable even though pseudopotentials.quantum-espresso.org is
  blocked by the network policy.
  QE inputs are written too, but QE needs UPF pseudopotentials which cannot be
  fetched here -- run those on a machine that can reach the QE pseudo library.
  The required files are listed in each generated input as a comment.

USAGE
  python 38_build_slabs.py                 # build everything into ../slabs/
  python 38_build_slabs.py --a 14.8        # override the HATCN 2D lattice const
"""
import os, argparse, shutil
import numpy as np
from ase import Atoms
from ase.build import bulk, surface
from ase.io import write as ase_write

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
STR = os.path.join(BASE, "structures")
OUT = os.path.join(BASE, "slabs")

def _find_cp2k_data():
    """CP2K ships GTH_POTENTIALS, BASIS_MOLOPT and dftd3.dat in its own data dir.
    Fall back to the copies inside the pyscf wheel (same format) if CP2K is not
    installed on the machine generating the inputs."""
    for root in ("/root/miniforge3/envs/cp2k/share/cp2k/data",
                 os.environ.get("CP2K_DATA_DIR", "")):
        if root and os.path.isfile(os.path.join(root, "GTH_POTENTIALS")):
            return root
    try:
        from pyscf.pbc.gto import pseudo as _ps
        return os.path.dirname(_ps.__file__)      # basis lives in a sibling dir
    except Exception:
        return None


CP2K_DATA = _find_cp2k_data()
PSEUDO_FILE = os.path.join(CP2K_DATA, "GTH_POTENTIALS") if CP2K_DATA else None
BASIS_FILE = os.path.join(CP2K_DATA, "BASIS_MOLOPT") if CP2K_DATA else None
D3_FILE = os.path.join(CP2K_DATA, "dftd3.dat") if CP2K_DATA else None

# MOLOPT basis + GTH-PBE potential per element. q-value must match the potential
# actually present in GTH_POTENTIALS; Ag q11 keeps 4d/5s, q19 adds the 4s4p
# semicore (safer for a metal, ~2x cost).
KINDS = {
    "H":  ("DZVP-MOLOPT-GTH", "GTH-PBE-q1"),
    "C":  ("DZVP-MOLOPT-GTH", "GTH-PBE-q4"),
    "N":  ("DZVP-MOLOPT-GTH", "GTH-PBE-q5"),
    "O":  ("DZVP-MOLOPT-GTH", "GTH-PBE-q6"),
    "F":  ("DZVP-MOLOPT-GTH", "GTH-PBE-q7"),
    "Li": ("DZVP-MOLOPT-SR-GTH", "GTH-PBE-q3"),
    "Al": ("DZVP-MOLOPT-SR-GTH", "GTH-PBE-q3"),
    "S":  ("DZVP-MOLOPT-GTH", "GTH-PBE-q6"),
    "Mo": ("DZVP-MOLOPT-SR-GTH", "GTH-PBE-q14"),
    "Ag": ("DZVP-MOLOPT-SR-GTH", "GTH-PBE-q11"),
}
# UPF names for the QE inputs (pslibrary PAW, PBE)
UPF = {
    "H": "H.pbe-kjpaw_psl.1.0.0.UPF",   "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
    "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF", "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "F": "F.pbe-n-kjpaw_psl.1.0.0.UPF", "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
    "Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF", "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Mo": "Mo.pbe-spn-kjpaw_psl.1.0.0.UPF", "Ag": "Ag.pbe-n-kjpaw_psl.1.0.0.UPF",
}


# ---------------------------------------------------------------- structures
def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0])
    sym = [l.split()[0] for l in L[2:2 + n]]
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in L[2:2 + n]])
    return sym, xyz


def hatcn_monolayer(a=15.0, vacuum=22.0):
    """Flat-lying HATCN in a hexagonal 2D cell.

    a is the in-plane lattice constant. HATCN spans 12.0 A N-to-N, so a must
    leave a van-der-Waals gap between the nitrile rims of neighbouring
    molecules: a = 15.0 A gives ~3.0 A N...N, about right for a physisorbed
    monolayer. THE VALUE IS AN ASSUMPTION -- job `hatcn_cellscan` scans it.
    """
    sym, x = read_xyz(os.path.join(STR, "HATCN.xyz"))
    x = x - x.mean(axis=0)
    # rotate the molecular plane onto xy (the molecule is planar to <1e-3 A)
    u, s, vt = np.linalg.svd(x - x.mean(axis=0))
    x = x @ vt.T                       # principal axes; 3rd is the normal
    x[:, 2] -= x[:, 2].mean()
    cell = np.array([[a, 0, 0],
                     [-a / 2, a * np.sqrt(3) / 2, 0],
                     [0, 0, vacuum]])
    x[:, 2] += vacuum / 2              # centre in the vacuum slab
    at = Atoms(symbols=sym, positions=x, cell=cell, pbc=[True, True, True])
    return at


def ag_on_hatcn(at, site="nitrile", height=2.3):
    """Place one Ag adatom.

    nitrile : along the C-N axis of a terminal nitrile, the site that gave
              E_b = 1.03 eV in the cluster model.
    ring    : above the centroid of the aromatic core (the inactive site, ~0 eV).
    bridge  : midway between two neighbouring nitrile N of the SAME molecule --
              the transition state candidate for surface diffusion.
    """
    at = at.copy()
    sym = at.get_chemical_symbols()
    p = at.get_positions()
    zc = p[:, 2].mean()
    ns = [i for i, s in enumerate(sym) if s == "N"]
    cen = p[:, :2].mean(axis=0)

    if site == "ring":
        pos = np.array([cen[0], cen[1], zc + 3.0])
    elif site == "nitrile":
        # outermost N = nitrile nitrogen; go outward in-plane and up
        i = max(ns, key=lambda k: np.linalg.norm(p[k, :2] - cen))
        d = p[i, :2] - cen
        d = d / np.linalg.norm(d)
        pos = np.array([p[i, 0] + height * d[0], p[i, 1] + height * d[1], zc + 1.2])
    elif site == "bridge":
        outer = sorted(ns, key=lambda k: -np.linalg.norm(p[k, :2] - cen))[:6]
        outer.sort(key=lambda k: np.arctan2(*(p[k, :2] - cen)[::-1]))
        m = 0.5 * (p[outer[0], :2] + p[outer[1], :2])
        d = m - cen
        d = d / np.linalg.norm(d)
        pos = np.array([m[0] + 1.2 * d[0], m[1] + 1.2 * d[1], zc + 1.5])
    else:
        raise ValueError(site)
    at.append(Atoms("Ag", positions=[pos])[0])
    return at


def lif_slab(a=4.026, layers=4, vacuum=18.0):
    b = bulk("LiF", "rocksalt", a=a)
    s = surface(b, (0, 0, 1), layers, vacuum=vacuum / 2)
    s = s.repeat((2, 2, 1))
    s.pbc = [True, True, True]
    return s


def ag_on_lif(slab, height=2.5):
    s = slab.copy()
    p = s.get_positions()
    sym = s.get_chemical_symbols()
    top_z = p[:, 2].max()
    tops = [i for i in range(len(s)) if p[i, 2] > top_z - 0.5 and sym[i] == "F"]
    i = tops[len(tops) // 2]
    s.append(Atoms("Ag", positions=[[p[i, 0], p[i, 1], p[i, 2] + height]])[0])
    return s


def neb_images(initial, final, n_images=7):
    """Linear interpolation of the Ag position only (substrate frozen).
    CI-NEB will relax them; this is just the starting band."""
    imgs = []
    p0 = initial.get_positions()[-1]
    p1 = final.get_positions()[-1]
    for k in range(n_images):
        at = initial.copy()
        pos = at.get_positions()
        pos[-1] = p0 + (p1 - p0) * k / (n_images - 1)
        at.set_positions(pos)
        imgs.append(at)
    return imgs


# ------------------------------------------------------------------- writers
CP2K_TMPL = """! {title}
! generated by scripts/38_build_slabs.py
! run:  cp2k.psmp -i {name}.inp -o {name}.out
&GLOBAL
  PROJECT {name}
  RUN_TYPE {run_type}
  PRINT_LEVEL MEDIUM
  WALLTIME {walltime}
&END GLOBAL
{motion}&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME {basis}
    POTENTIAL_FILE_NAME {pseudo}
    &QS
      EPS_DEFAULT 1.0E-12
      METHOD GPW
    &END QS
    &MGRID
      ! CONVERGE THESE FIRST -- see jobs/cutoff_scan
      CUTOFF {cutoff}
      REL_CUTOFF {rel_cutoff}
      NGRIDS 5
    &END MGRID
    &POISSON
      PERIODIC XY
      POISSON_SOLVER {poisson}
    &END POISSON
{scf}    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
{vdw}    &END XC
{uks}  &END DFT
  &SUBSYS
    &CELL
      A {A}
      B {B}
      C {C}
      PERIODIC XYZ
    &END CELL
    &COORD
{coord}    &END COORD
{kinds}  &END SUBSYS
&END FORCE_EVAL
"""

# Ag present -> partial occupations are real; smearing + diagonalisation.
SCF_METAL = """    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 60
      ADDED_MOS {added_mos}
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE 300
      &END SMEAR
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.2
        NBUFFER 8
      &END MIXING
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
"""

# Gapped, closed-shell system -> orbital transformation, no smearing.
# Applying Fermi smearing to an insulator (as a metal template would) makes the
# SCF oscillate at the 1e-6 level instead of converging: measured on lif001,
# which needed 27+ noisy iterations and 609 s on 4 cores that way.
SCF_GAP = """    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 30
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
"""

VDW_D3 = """      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
"""

MOTION_GEOOPT = """&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 200
    MAX_FORCE 4.5E-4
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      ! freeze the substrate, relax the adatom only, for a first pass
      LIST {frozen}
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION
"""


def write_cp2k(at, name, out_dir, title="", run_type="ENERGY_FORCE",
               cutoff=500, rel_cutoff=60, uks=False, d3=True,
               frozen=None, walltime=86000):
    os.makedirs(out_dir, exist_ok=True)
    cell = at.get_cell()
    coord = "".join(f"      {s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}\n"
                    for s, p in zip(at.get_chemical_symbols(), at.get_positions()))
    kinds = ""
    for el in sorted(set(at.get_chemical_symbols())):
        b, q = KINDS[el]
        kinds += (f"    &KIND {el}\n      BASIS_SET {b}\n"
                  f"      POTENTIAL {q}\n    &END KIND\n")
    n_ag = sum(1 for s in at.get_chemical_symbols() if s == "Ag")
    motion = ""
    if run_type == "GEO_OPT" and frozen:
        motion = MOTION_GEOOPT.format(frozen=frozen)
    txt = CP2K_TMPL.format(
        title=title or name, name=name, run_type=run_type, walltime=walltime,
        basis=os.path.basename(BASIS_FILE) if BASIS_FILE else "BASIS_MOLOPT",
        pseudo=os.path.basename(PSEUDO_FILE) if PSEUDO_FILE else "GTH_POTENTIALS",
        cutoff=cutoff, rel_cutoff=rel_cutoff,
        poisson="ANALYTIC" if at.pbc[2] else "MT",
        scf=(SCF_METAL.format(added_mos=max(40, 20 * n_ag)) if n_ag else SCF_GAP),
        vdw=VDW_D3 if d3 else "",
        uks="    UKS .TRUE.\n    MULTIPLICITY 2\n" if uks else "",
        A=" ".join(f"{v:12.6f}" for v in cell[0]),
        B=" ".join(f"{v:12.6f}" for v in cell[1]),
        C=" ".join(f"{v:12.6f}" for v in cell[2]),
        coord=coord, kinds=kinds, motion=motion)
    open(os.path.join(out_dir, f"{name}.inp"), "w").write(txt)


QE_TMPL = """! {title}
! generated by scripts/38_build_slabs.py
! NOTE: needs UPF pseudopotentials in pseudo_dir. Fetch from
!   https://pseudopotentials.quantum-espresso.org/upf_files/<name>
! required here: {upflist}
&CONTROL
  calculation = '{calc}'
  prefix = '{name}'
  outdir = './out'
  pseudo_dir = './pseudo'
  tprnfor = .true.
  tstress = .false.
/
&SYSTEM
  ibrav = 0
  nat = {nat}
  ntyp = {ntyp}
  ecutwfc = {ecutwfc}
  ecutrho = {ecutrho}
  occupations = 'smearing'
  smearing = 'mv'
  degauss = 0.01
  vdw_corr = 'grimme-d3'
  assume_isolated = '2D'
/
&ELECTRONS
  conv_thr = 1.0d-6
  mixing_beta = 0.3
  electron_maxstep = 200
/
&IONS
  ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
{species}
CELL_PARAMETERS angstrom
{cellblock}
ATOMIC_POSITIONS angstrom
{positions}
K_POINTS automatic
{kmesh} 1 0 0 0
"""


def write_qe(at, name, out_dir, title="", calc="scf",
             ecutwfc=60, ecutrho=480, kmesh="2 2"):
    os.makedirs(out_dir, exist_ok=True)
    els = sorted(set(at.get_chemical_symbols()))
    from ase.data import atomic_masses, atomic_numbers
    species = "".join(f"  {e} {atomic_masses[atomic_numbers[e]]:.4f} {UPF[e]}\n"
                      for e in els)
    cellblock = "".join("  " + " ".join(f"{v:12.6f}" for v in r) + "\n"
                        for r in at.get_cell())
    positions = "".join(f"  {s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n"
                        for s, p in zip(at.get_chemical_symbols(), at.get_positions()))
    txt = QE_TMPL.format(
        title=title or name, name=name, calc=calc, nat=len(at), ntyp=len(els),
        ecutwfc=ecutwfc, ecutrho=ecutrho, species=species, cellblock=cellblock,
        positions=positions, kmesh=kmesh,
        upflist=", ".join(UPF[e] for e in els))
    open(os.path.join(out_dir, f"{name}.in"), "w").write(txt)


# ---------------------------------------------------------------------- main
def main(a_hatcn=15.0):
    os.makedirs(OUT, exist_ok=True)
    # Inputs reference the data files by bare name; CP2K resolves them against
    # $CP2K_DATA_DIR. Write a run helper that sets it rather than copying ~10 MB
    # of tables into the repo.
    if CP2K_DATA:
        open(os.path.join(OUT, "run_cp2k.sh"), "w").write(
            "#!/bin/bash\n"
            "# usage: ./run_cp2k.sh <jobdir>/<name>.inp\n"
            "source /root/miniforge3/etc/profile.d/conda.sh\n"
            "conda activate cp2k\n"
            f"export CP2K_DATA_DIR={CP2K_DATA}\n"
            'export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}\n'
            'inp="$1"; cd "$(dirname "$inp")"\n'
            'cp2k.psmp -i "$(basename "$inp")" -o "$(basename "${inp%.inp}").out"\n')
        os.chmod(os.path.join(OUT, "run_cp2k.sh"), 0o755)

    built = []

    ml = hatcn_monolayer(a=a_hatcn)
    n_sub = len(ml)
    frozen_sub = f"1..{n_sub}"

    jobs = [
        ("hatcn_ml", ml, "ENERGY_FORCE", False, None,
         "HATCN monolayer, clean reference"),
        ("hatcn_ml_ag", ag_on_hatcn(ml, "nitrile"), "GEO_OPT", True, frozen_sub,
         "Ag adatom on nitrile N (cluster value 1.03 eV)"),
        ("hatcn_ml_ag_ring", ag_on_hatcn(ml, "ring"), "GEO_OPT", True, frozen_sub,
         "Ag adatom over the aromatic core (cluster value ~0 eV)"),
        ("hatcn_ml_ag_bridge", ag_on_hatcn(ml, "bridge"), "GEO_OPT", True, frozen_sub,
         "Ag at the bridge site -- diffusion transition-state candidate"),
    ]
    lif = lif_slab()
    jobs += [
        ("lif001", lif, "ENERGY_FORCE", False, None,
         "LiF(001) 2x2 4-layer slab, clean"),
        ("lif001_ag", ag_on_lif(lif), "GEO_OPT", True, f"1..{len(lif)}",
         "Ag adatom on LiF(001) F-top (cluster value 0.25 eV)"),
    ]

    for name, at, rt, uks, frozen, desc in jobs:
        d = os.path.join(OUT, name)
        os.makedirs(d, exist_ok=True)
        ase_write(os.path.join(d, f"{name}.xyz"), at)
        ase_write(os.path.join(d, f"{name}.cif"), at)
        write_cp2k(at, name, d, title=desc, run_type=rt, uks=uks, frozen=frozen)
        write_qe(at, name, d, title=desc,
                 calc="relax" if rt == "GEO_OPT" else "scf",
                 kmesh="1 1" if "hatcn" in name else "3 3")
        built.append((name, len(at), desc))

    # NEB band: nitrile -> bridge -> next nitrile
    ini = ag_on_hatcn(ml, "nitrile")
    fin = ag_on_hatcn(ml, "bridge")
    d = os.path.join(OUT, "hatcn_neb")
    os.makedirs(d, exist_ok=True)
    imgs = neb_images(ini, fin, 7)
    for i, im in enumerate(imgs):
        ase_write(os.path.join(d, f"image_{i:02d}.xyz"), im)
    ase_write(os.path.join(d, "band.xyz"), imgs)
    built.append(("hatcn_neb", len(ini), "7-image CI-NEB starting band"))

    # convergence-test inputs: cutoff scan on the cheap system
    d = os.path.join(OUT, "cutoff_scan")
    for c in (300, 400, 500, 600, 700):
        write_cp2k(lif, f"lif_cut{c}", d, run_type="ENERGY_FORCE",
                   cutoff=c, title=f"LiF(001) cutoff convergence, {c} Ry")
    built.append(("cutoff_scan", len(lif), "MGRID CUTOFF 300-700 Ry on LiF(001)"))

    print(f"{'job':<22}{'atoms':>7}  description")
    for n, k, dsc in built:
        print(f"{n:<22}{k:>7}  {dsc}")
    print(f"\nwritten under {OUT}/")
    if CP2K_DATA:
        print(f"CP2K data dir: {CP2K_DATA}")
        print(f"run a job with:  {OUT}/run_cp2k.sh <job>/<name>.inp")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=float, default=15.0,
                    help="HATCN monolayer 2D lattice constant (A)")
    args = ap.parse_args()
    main(args.a)
