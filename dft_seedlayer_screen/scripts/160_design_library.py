"""Design library for a charge-generation-side Ag seed: HAT and phenanthroline scaffolds.

    python scripts/160_design_library.py  ->  runs/design/library.json

Role fixed first (notes/DESIGN_SPEC_CGL_SEED.md): the seed sits on an
n-doped ETL under a thin Ag cathode, as the p-side of the charge-generation
junction (n-ETL / seed / Ag), like HATCN. So the scaffold must stay a deep
acceptor; what may change is the periphery.

Only molecules with a known route are generated:

HAT core (hexaazatriphenylene), three pyrazine "units" around a benzene:
  route A  hexaketocyclohexane + 1,2-diamines   (HATCN: diaminomaleonitrile)
           units: CN (C#N,C#N), benzo family (o-phenylenediamines: H4, F2,
           Cl2, (CN)2, F4, mono-CN), pyrido[3,4-b] (3,4-diaminopyridine)
  route B  hexaaminobenzene + 1,2-diketones
           units: H (glyoxal), Me (diacetyl), Ph (benzil), C6F5
           (decafluorobenzil), 4-CN-Ph (4,4'-dicyanobenzil), 2-py
           (2,2'-pyridil), phen (1,10-phenanthroline-5,6-dione)
  Three units with repetition. A molecule mixing route-A and route-B units
  has no one-pot route and is flagged; two different units of one route
  give a statistical mixture that must be separated (flagged, not dropped).

Phenanthroline:
  dppz family  phen-5,6-dione + the route-A diamines (keeps the phen chelate,
               adds an acceptor phenazine)
  Bphen family 4,7-diaryl-1,10-phenanthrolines (Suzuki on 4,7-dichloro-phen)

References (role anchors, computed the same way): HATCN, F4TCNQ, F6-TCNNQ,
Bphen, BCP, TPBi.
"""
import itertools
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, RDConfig, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "runs", "design", "library.json")

# a unit is the pyrazine C-C pair, written between the two pyrazine N atoms;
# {a},{b},{c} are ring-closure numbers unique to the unit instance
PAIR = {  # route, substituent smiles (R)
    "CN":   ("A", "C#N"),
    "H":    ("B", None),
    "Me":   ("B", "C"),
    "Ph":   ("B", "c{a}ccccc{a}"),
    "C6F5": ("B", "c{a}c(F)c(F)c(F)c(F)c{a}F"),
    "CNPh": ("B", "c{a}ccc(C#N)cc{a}"),
    "2py":  ("B", "c{a}ccccn{a}"),
}
FUSED = {  # route, fused ring; first and last atoms are the pyrazine carbons
    "benzo":   ("A", "c{a}ccccc{a}"),
    "F2benzo": ("A", "c{a}cc(F)c(F)cc{a}"),
    "Cl2benzo": ("A", "c{a}cc(Cl)c(Cl)cc{a}"),
    "CN2benzo": ("A", "c{a}cc(C#N)c(C#N)cc{a}"),
    "F4benzo": ("A", "c{a}c(F)c(F)c(F)c(F)c{a}"),
    "CNbenzo": ("A", "c{a}ccc(C#N)cc{a}"),
    "pyr34":   ("A", "c{a}ccncc{a}"),
    "phen":    ("B", "c{a}c{b}cccnc{b}c{b}ncccc{b}c{a}"),
}
UNITS = {**{k: ("pair",) + v for k, v in PAIR.items()},
         **{k: ("fused",) + v for k, v in FUSED.items()}}


class Rings:
    def __init__(self, start=10):
        self.n = start

    def take(self):
        self.n += 1
        return f"%{self.n}"


def unit_smiles(name, rings):
    kind, route, s = UNITS[name]
    if kind == "pair":
        if s is None:
            return "cc"
        r1 = s.format(a=rings.take())
        r2 = s.format(a=rings.take())
        return f"c({r1})c({r2})"
    return s.format(a=rings.take(), b=rings.take())


def hat(u1, u2, u3):
    rings = Rings()
    a, b, c = (unit_smiles(u, rings) for u in (u1, u2, u3))
    # central benzene A..F; pyrazines (A,B), (C,D), (E,F)
    smi = f"c12c3n{a}nc3c3n{b}nc3c1n{c}n2"
    return smi


def dppz(u):
    rings = Rings()
    a = unit_smiles(u, rings)
    # phen with the pyrazine fused at C5-C6
    return f"c12c3cccnc3c3ncccc3c1n{a}n2"


ARYL = {"Ph": "c1ccccc1", "C6F5": "c1c(F)c(F)c(F)c(F)c1F", "CNPh": "c1ccc(C#N)cc1",
        "CF3Ph": "c1ccc(C(F)(F)F)cc1", "3py": "c1cccnc1", "4py": "c1ccncc1",
        "Mes": "c1c(C)cc(C)cc1C", "Naph": "c1cccc2ccccc12", "35CF3Ph": "c1cc(C(F)(F)F)cc(C(F)(F)F)c1"}


def bphen(ar):
    s = ARYL[ar]
    s1 = s.replace("1", "%21")
    s2 = s.replace("1", "%22")
    return f"c1cnc2c(c1{s1})ccc1c(ccnc12){s2}"


REFERENCES = {
    "HATCN": "N#Cc1nc2c3nc(C#N)c(C#N)nc3c3nc(C#N)c(C#N)nc3c2nc1C#N",
    "F4TCNQ": "N#CC(C#N)=C1C(F)=C(F)C(=C(C#N)C#N)C(F)=C1F",
    "F6-TCNNQ": "N#CC(C#N)=C1C(F)=C(F)C2=C(C(F)=C(F)C(=C(C#N)C#N)C2=C1)",
    "Bphen": "c1ccc(cc1)-c1ccnc2c1ccc1c(ccnc12)-c1ccccc1",
    "BCP": "Cc1cc(-c2ccccc2)c2ccc3c(-c4ccccc4)cc(C)nc3c2n1",
    "TPBi": "c1ccc(cc1)-n1c(nc2ccccc12)-c1cc(cc(c1)-c1nc2ccccc2n1-c1ccccc1)-c1nc2ccccc2n1-c1ccccc1",
}


def describe(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
    heavy = m.GetNumHeavyAtoms()
    return {"smiles": Chem.MolToSmiles(m), "formula": rdMolDescriptors.CalcMolFormula(m),
            "MW": round(Descriptors.MolWt(m), 1), "heavy": heavy,
            "HBD": rdMolDescriptors.CalcNumHBD(m),
            "nrot": rdMolDescriptors.CalcNumRotatableBonds(m),
            "sym_frac": round(len(set(ranks)) / heavy, 3),
            "SA": round(sascorer.calculateScore(m), 2),
            "n_CN": smi.count("C#N"), "n_aromN": sum(1 for a in m.GetAtoms()
                                                       if a.GetSymbol() == "N" and a.GetIsAromatic())}


def main():
    lib = {}
    for name, smi in REFERENCES.items():
        d = describe(smi)
        d.update(family="reference", route="known", units=[name])
        lib[name] = d
    names = sorted(UNITS)
    for combo in itertools.combinations_with_replacement(names, 3):
        smi = hat(*combo)
        d = describe(smi)
        if d is None:
            print("bad", combo)
            continue
        routes = {UNITS[u][1] for u in combo}
        distinct = len(set(combo))
        d.update(family="HAT", units=list(combo),
                 route="mixed A+B (no one-pot route)" if len(routes) > 1 else
                 f"{routes.pop()}, {'homoleptic' if distinct == 1 else 'statistical mixture'}",
                 synth_class=3 if len(routes) > 1 else (1 if distinct == 1 else 2))
        lib["HAT-" + "-".join(combo)] = d
    for u in names:
        if UNITS[u][1] != "A" and u not in ("CN",):
            continue
        d = describe(dppz(u))
        if d:
            d.update(family="dppz", units=[u], route="phen-5,6-dione + diamine", synth_class=1)
            lib["dppz-" + u] = d
    for ar in ARYL:
        d = describe(bphen(ar))
        if d:
            d.update(family="Bphen", units=[ar], route="4,7-dichloro-phen + ArB(OH)2", synth_class=1)
            lib["Bphen-" + ar] = d
    # dedupe by canonical smiles
    seen, out = {}, {}
    for k, v in lib.items():
        if v["smiles"] in seen:
            continue
        seen[v["smiles"]] = k
        out[k] = v
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1)
    fam = {}
    for v in out.values():
        fam[v["family"]] = fam.get(v["family"], 0) + 1
    print(len(out), fam)
    print("HATCN check:", out["HAT-CN-CN-CN"]["smiles"] == out["HATCN"]["smiles"]
          if "HAT-CN-CN-CN" in out else "deduped into reference")


if __name__ == "__main__":
    main()
