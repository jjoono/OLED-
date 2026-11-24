"""First look at non-HAT scaffolds built on HATCN's design rules.

    python scripts/176_new_scaffolds.py   ->  runs/design/new_scaffolds.json

What makes HATCN the seed (scripts 150, 166, 174): a strong acceptor (E_b far
above E* = 0.30 eV), the highest anchor density of any material screened
(12 lone-pair N on a 317 A^2 molecule, D = 37.8), and a wide optical gap. Its
weakness is morphology: a flat, rigid, D3h disc that crystallises.
Each idea below keeps the anchors dense and the core electron-poor and changes
one thing. Properties are the cheap tier of the design pipeline (161: GFN2
geometry, IPEA1-xTB IP/EA, GFN2 gap, crystallisation-risk score) plus D from
script 174 assuming every anchor traps (true for strong acceptors; checked
later at PBE0). EA and gap are read relative to HATCN only.
"""
import importlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
OUT = os.path.join(ROOT, "runs", "design", "new_scaffolds.json")

IDEAS = {
    # name: (SMILES, idea)
    "HATCN": ("N#Cc1nc2c3nc(C#N)c(C#N)nc3c3nc(C#N)c(C#N)nc3c2nc1C#N", "reference"),
    "TCNPz": ("N#Cc1nc(C#N)c(C#N)nc1C#N", "smallest unit: tetracyanopyrazine (too volatile alone)"),
    "TCN-Qx": ("N#Cc1nc2cc(C#N)c(C#N)cc2nc1C#N", "two rings, quinoxaline tetracarbonitrile"),
    "TCN-PQ": ("N#Cc1nc2cc3nc(C#N)c(C#N)nc3cc2nc1C#N", "linear three rings, lower symmetry than HAT"),
    "CN3-heptazine": ("N#Cc1nc2nc(C#N)nc3nc(C#N)nc(n1)n23", "tri-s-triazine core: most N per area"),
    "triptycene-HAT": ("N#Cc1nc2cc3c(cc2nc1C#N)C1c2cc4nc(C#N)c(C#N)nc4cc2C3c2cc3nc(C#N)c(C#N)nc3cc21",
                       "3D: three dicyanopyrazino arms on a triptycene, cannot stack like a disc"),
    "bi-TCNQx": ("N#Cc1nc2cc(-c3cc4nc(C#N)c(C#N)nc4cc3C#N)c(C#N)cc2nc1C#N",
                 "twisted biaryl dimer of TCN-Qx: dense anchors, non-planar"),
    # heptazine dimers (2026-11-24): fix the monomer's MW 248 (volatile) and
    # flat-disc packing while keeping its anchor density
    "biHz": ("N#Cc1nc2nc(-c4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc3nc(C#N)nc(n1)n23",
             "heptazine dimer, direct C-C bond"),
    "Hz-NCN-Hz": ("N#Cc1nc2nc(N(C#N)c4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc3nc(C#N)nc(n1)n23",
                  "cyanamide-bridged dimer (melon-type N bridge, bent, extra nitrile anchor)"),
    "Hz-O-Hz": ("N#Cc1nc2nc(Oc4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc3nc(C#N)nc(n1)n23",
                "ether-bridged dimer (bent, flexible)"),
    "Hz-C6F4-Hz": ("N#Cc1nc2nc(-c7c(F)c(F)c(-c4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)c(F)c7F)nc3nc(C#N)nc(n1)n23",
                   "tetrafluorophenylene-bridged dimer (F forces a twist)"),
    # second round (2026-11-24): bridges and oligomers around Hz-NCN-Hz
    "Hz-C(CN)2-Hz": ("N#Cc1nc2nc(C(C#N)(C#N)c4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc3nc(C#N)nc(n1)n23",
                     "sp3 dicyanomethylene bridge: tetrahedral kink, two extra nitrile anchors"),
    "Hz-SO2-Hz": ("N#Cc1nc2nc(S(=O)(=O)c4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc3nc(C#N)nc(n1)n23", "sulfone bridge: strongly bent, electron-withdrawing"),
    "Hz-S-Hz": ("N#Cc1nc2nc(Sc4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc3nc(C#N)nc(n1)n23", "thioether bridge: bent, soft S"),
    "Hz3N": ("N(c1nc2nc(C#N)nc3nc(C#N)nc(n1)n23)(c1nc2nc(C#N)nc3nc(C#N)nc(n1)n23)c1nc2nc(C#N)nc3nc(C#N)nc(n1)n23",
             "tri-heptazinyl amine: C3 propeller, MW ~700"),
    "Hz-NCN-Hz-NCN-Hz": ("N#Cc1nc2nc(N(C#N)c4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc3nc(N(C#N)c4nc5nc(C#N)nc6nc(C#N)nc(n4)n56)nc(n1)n23",
                         "linear trimer with two cyanamide bridges"),
}


def main():
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    s161 = importlib.import_module("161_design_props")
    t = importlib.import_module("174_two_tier_descriptor")
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for name, (smi, idea) in IDEAS.items():
        if name in out:
            continue
        m = Chem.MolFromSmiles(smi)
        heavy = m.GetNumHeavyAtoms()
        ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
        d = {"smiles": smi, "MW": round(Descriptors.MolWt(m), 1),
             "nrot": rdMolDescriptors.CalcNumRotatableBonds(m),
             "sym_frac": round(len(set(ranks)) / heavy, 3)}
        p = s161.work((name, d))
        n, sa = t.anchors_sasa(t.mol_from_smiles(smi))
        p.update(idea=idea, formula=rdMolDescriptors.CalcMolFormula(m), MW=d["MW"],
                 n_anchor=n, SASA=round(sa, 1), D_if_all_trap=round(1e3 * n / sa, 1))
        out[name] = p
        json.dump(out, open(OUT, "w"), indent=1)
        print(f"{name:<15} {p['formula']:<10} MW {d['MW']:6.1f}  EA {p.get('EA')}  gap {p.get('gap_gfn2')}  "
              f"risk {p.get('cryst_risk')}  3D {p.get('threeD')}  anchors {n}  D {p['D_if_all_trap']}", flush=True)


if __name__ == "__main__":
    main()
