"""Ag-M dimer binding energies (M = Al, Mg, Cu, Ag) at PBE-D3BJ/def2-SVP with CP,
same level as the seed-layer screening -> directly comparable to Ag2 = 1.86 eV.
"""
import psi4, os, json

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")
psi4.set_memory("4 GB")
psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_dimers.out"), False)
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "uks",
                  "maxiter": 200, "guess": "sad"})
M = "pbe-d3bj"
H2EV = 27.211386

# (partner, dimer multiplicity, start distance A)
cases = {"Al": (2, 2.60), "Mg": (2, 2.90), "Cu": (2, 2.40), "Ag": (1, 2.60)}
# dimer mult: AgAl -> singlet(1) or triplet? AgAl ground state is 1Sigma+ -> mult 1
mults = {"Al": 1, "Mg": 2, "Cu": 1, "Ag": 1}
atom_mult = {"Al": 2, "Mg": 1, "Cu": 2, "Ag": 2}

res = {}
for m, (amult, r0) in cases.items():
    dm = mults[m]
    mol = psi4.geometry(f"0 {dm}\nAg 0 0 0\n{m} 0 0 {r0}\nsymmetry c1\n")
    e_d = psi4.optimize(M, molecule=mol)
    r = abs(mol.geometry().np[1][2] - mol.geometry().np[0][2]) * 0.52917721067
    # CP fragments at dimer geometry
    g = mol.geometry().np * 0.52917721067
    s = f"0 2\nAg {g[0][0]} {g[0][1]} {g[0][2]}\nGh({m}) {g[1][0]} {g[1][1]} {g[1][2]}\nsymmetry c1\n"
    e_ag = psi4.energy(M, molecule=psi4.geometry(s))
    s = f"0 {atom_mult[m]}\nGh(Ag) {g[0][0]} {g[0][1]} {g[0][2]}\n{m} {g[1][0]} {g[1][1]} {g[1][2]}\nsymmetry c1\n"
    e_m = psi4.energy(M, molecule=psi4.geometry(s))
    be = (e_ag + e_m - e_d) * H2EV
    res[f"Ag{m}"] = {"BE_eV": be, "r_A": r}
    print(f"Ag-{m}: BE = {be:.3f} eV, r = {r:.3f} A", flush=True)
    psi4.core.clean()

json.dump(res, open(os.path.join(RUNS, "dimer_BE.json"), "w"), indent=2)
