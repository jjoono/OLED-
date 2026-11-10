"""GFN2 -> PBE0 map for the Ag binding energy, with uncertainty (Delta-ML).

    python scripts/163_delta_ml.py           # fit, LOO, write runs/design/dml_model.json
    python scripts/163_delta_ml.py --apply   # + predict the library (runs/design/sites.jsonl)

Why a map is needed: GFN2 finds the right site but overbinds anything with
nitrogen -- PhCN 0.95 vs PBE0 0.23 eV, HATCN 3.29 vs 1.63. The overbinding is
not a constant factor: it is largest on weak acceptors (a nitrile N that PBE0
barely binds) and smaller on deep acceptors, where the electron transfer GFN2
sees is real. So the map is fed what separates those two cases:

  x1  GFN2 E_b at the deepest site
  x2  GFN2 q(Ag1)                       (how much charge actually moved)
  x3  IPEA1-xTB vertical EA             (could the molecule take it)

Training set: every molecule with both a PBE0 multi-site (or survey-confirmed)
E_b and a GFN2 survey in this project -- 22 common OLED materials (script 150)
and 13 small probes (script 141). EA computed here the same way as script 161.

Two models, compared by leave-one-out:
  ridge  linear in (x1, x2, x3, x1*x3)
  GP     ARD-RBF + white noise on standardised x, hyperparameters by marginal
         likelihood; gives the predictive sigma used as the uncertainty.
The predictive sigma is the honest part: the library's best HAT derivatives
have EA above every training point but three (HATCN, F4TCNQ, F6-TCNNQ), so
their predictions are extrapolations, and sigma says so.
"""
import json
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
from scipy.optimize import minimize  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.join(HERE, "..")
D = os.path.join(ROOT, "runs", "design")
TRAIN = os.path.join(D, "dml_train.json")
MODEL = os.path.join(D, "dml_model.json")
PRED = os.path.join(D, "dml_pred.json")

SMILES = {  # training molecules without a GFN2 structure in structures/common
    "p-bPPhenB": "c1ccc(cc1)-c1cc(-c2ccc(cc2)-c2cc(-c3ccccc3)nc3c2ccc2cccnc23)c2ccc3cccnc3c2n1",
    "B3PyMPM": "Cc1nc(-c2cc(-c3cccnc3)cc(-c3cccnc3)c2)cc(-c2cc(-c3cccnc3)cc(-c3cccnc3)c2)n1",
    "Liq": "[Li]Oc1cccc2cccnc12",
    "pyridine": "c1ccncc1", "triazine": "c1ncncn1", "BTD": "c1ccc2nsnc2c1",
    "Me3P=O": "CP(C)(C)=O", "Ph3P=O": "O=P(c1ccccc1)(c1ccccc1)c1ccccc1",
    "thiophene": "c1ccsc1", "benzene": "c1ccccc1", "PhCN": "N#Cc1ccccc1",
    "DMABN": "CN(C)c1ccc(C#N)cc1", "pDCNB": "N#Cc1ccc(C#N)cc1", "oDCNB": "N#Cc1ccccc1C#N",
    "PhCz": "c1ccc(-n2c3ccccc3c3ccccc32)cc1", "TPA": "c1ccc(N(c2ccccc2)c2ccccc2)cc1",
}
SMALL = ["pyridine", "triazine", "BTD", "Me3P=O", "Ph3P=O", "thiophene", "benzene", "PhCN",
         "DMABN", "pDCNB", "oDCNB", "PhCz", "TPA"]


def ea_from_xyz(S, X):
    s161 = __import__("importlib").import_module("161_design_props")
    from ase import Atoms
    Z = np.array(Atoms(S).numbers)
    X = np.asarray(X, float)
    e0 = s161.sp("IPEA1-xTB", Z, X, 0, 0).get("energy")
    em = s161.sp("IPEA1-xTB", Z, X, -1, 1).get("energy")
    return float((e0 - em) * s161.HA - s161.SHIFT)


def ea_from_smiles(name, smi):
    s161 = __import__("importlib").import_module("161_design_props")
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    m = Chem.MolFromSmiles(smi)
    d = {"smiles": smi, "nrot": 0, "sym_frac": 1.0, "MW": Descriptors.MolWt(m)}
    r = s161.work((name, d))
    return r["EA"]


def build_train():
    from pathgeom import read_xyz
    props = {}
    p = os.path.join(D, "props.jsonl")
    if os.path.exists(p):
        props = {json.loads(l)["name"]: json.loads(l) for l in open(p) if l.strip()}
    common = json.load(open(os.path.join(ROOT, "runs", "common_survey_xtb.json")))
    old = json.load(open(os.path.join(ROOT, "runs", "xtb_survey.json")))
    rank = {r["name"]: r for r in json.load(open(os.path.join(ROOT, "runs", "common_ranking.json")))}
    status = {r["candidate"]: r for r in json.load(open(os.path.join(ROOT, "runs", "screening_status.json")))}
    rows = []
    for name in list(rank) + SMALL:
        if name in rank:
            pbe0 = rank[name]["E_b"]
        else:
            pbe0 = status[name]["E_b_deepest"]
        if name in common:
            g = common[name]
            eb_g, q1 = g["sites"][g["deepest_site"]]["E_b"], g["q_Ag1"]
        else:
            g = old[name]
            eb_g, q1 = g["E_b_deepest"], g["q_Ag1"]
        if name in props and "EA" in props[name]:
            ea = props[name]["EA"]
        elif os.path.exists(os.path.join(ROOT, "structures", "common", f"{name}.xyz")):
            S, X = read_xyz(os.path.join(ROOT, "structures", "common", f"{name}.xyz"))
            ea = ea_from_xyz(S, X)
        else:
            ea = ea_from_smiles(name, SMILES[name])
        rows.append({"name": name, "pbe0": round(float(pbe0), 4), "eb_gfn2": round(float(eb_g), 4),
                     "q1": round(float(q1), 4), "EA": round(float(ea), 3)})
        print(rows[-1], flush=True)
    json.dump(rows, open(TRAIN, "w"), indent=1)
    return rows


def feats(eb, q, ea, kind):
    eb, q, ea = map(np.atleast_1d, (eb, q, ea))
    if kind == "ridge":
        return np.column_stack([eb, q, ea, eb * ea])
    return np.column_stack([eb, q, ea])


class Ridge:
    def __init__(self, lam=1e-2):
        self.lam = lam

    def fit(self, X, y):
        self.mu, self.sd = X.mean(0), X.std(0) + 1e-9
        Z = np.column_stack([np.ones(len(X)), (X - self.mu) / self.sd])
        R = self.lam * np.eye(Z.shape[1])
        R[0, 0] = 0
        self.w = np.linalg.solve(Z.T @ Z + R, Z.T @ y)
        self.s = float(np.sqrt(np.mean((Z @ self.w - y) ** 2)))
        return self

    def predict(self, X):
        Z = np.column_stack([np.ones(len(X)), (X - self.mu) / self.sd])
        return Z @ self.w, np.full(len(X), self.s)


class GP:
    """Zero-mean GP on y - linear prior mean in x1 (so it falls back to a sane
    trend far from data rather than to the training mean)."""

    def k(self, A, B, th):
        ls = np.exp(th[:A.shape[1]])
        sf2 = np.exp(2 * th[A.shape[1]])
        d = ((A[:, None, :] - B[None, :, :]) / ls) ** 2
        return sf2 * np.exp(-0.5 * d.sum(-1))

    def nll(self, th, Z, r):
        sn2 = np.exp(2 * th[-1])
        K = self.k(Z, Z, th) + (sn2 + 1e-8) * np.eye(len(Z))
        try:
            L = np.linalg.cholesky(K)
        except np.linalg.LinAlgError:
            return 1e10
        a = np.linalg.solve(L.T, np.linalg.solve(L, r))
        return 0.5 * r @ a + np.log(np.diag(L)).sum()

    def fit(self, X, y):
        self.mu, self.sd = X.mean(0), X.std(0) + 1e-9
        Z = (X - self.mu) / self.sd
        self.lin = np.polyfit(Z[:, 0], y, 1)
        r = y - np.polyval(self.lin, Z[:, 0])
        best = None
        for ls0 in (0.0, 0.7):
            th0 = np.r_[np.full(Z.shape[1], ls0), np.log(r.std() + 1e-3), np.log(0.05)]
            o = minimize(self.nll, th0, args=(Z, r), method="L-BFGS-B",
                         bounds=[(-2, 3)] * Z.shape[1] + [(-5, 2), (np.log(0.01), np.log(0.5))])
            if best is None or o.fun < best.fun:
                best = o
        self.th, self.Z, self.r = best.x, Z, r
        sn2 = np.exp(2 * self.th[-1])
        self.K = self.k(Z, Z, self.th) + (sn2 + 1e-8) * np.eye(len(Z))
        self.a = np.linalg.solve(self.K, r)
        return self

    def predict(self, X):
        Z = (X - self.mu) / self.sd
        Ks = self.k(Z, self.Z, self.th)
        m = Ks @ self.a + np.polyval(self.lin, Z[:, 0])
        v = np.exp(2 * self.th[len(self.mu)]) - np.einsum("ij,ji->i", Ks, np.linalg.solve(self.K, Ks.T))
        return m, np.sqrt(np.maximum(v, 0) + np.exp(2 * self.th[-1]))


def loo(rows, kind):
    y = np.array([r["pbe0"] for r in rows])
    X = feats([r["eb_gfn2"] for r in rows], [r["q1"] for r in rows], [r["EA"] for r in rows], kind)
    pred, sig = np.zeros(len(y)), np.zeros(len(y))
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        mdl = (Ridge() if kind == "ridge" else GP()).fit(X[m], y[m])
        p, s = mdl.predict(X[i:i + 1])
        pred[i], sig[i] = p[0], s[0]
    return y, pred, sig


def main():
    rows = json.load(open(TRAIN)) if os.path.exists(TRAIN) and "--rebuild" not in sys.argv else build_train()
    out = {"n_train": len(rows)}
    yg = np.array([r["eb_gfn2"] for r in rows])
    y = np.array([r["pbe0"] for r in rows])
    out["raw_gfn2"] = {"spearman": round(float(spearmanr(yg, y)[0]), 3)}
    for kind in ("ridge", "gp"):
        y, p, s = loo(rows, kind)
        z = (y - p) / s
        out[kind] = {"loo_mae": round(float(np.abs(y - p).mean()), 3),
                     "loo_rmse": round(float(np.sqrt(((y - p) ** 2).mean())), 3),
                     "loo_spearman": round(float(spearmanr(y, p)[0]), 3),
                     "frac_within_1sigma": round(float((np.abs(z) < 1).mean()), 2),
                     "worst": sorted(((r["name"], round(float(yy - pp), 3)) for r, yy, pp in
                                      zip(rows, y, p)), key=lambda t: -abs(t[1]))[:5],
                     "loo": {r["name"]: [round(float(pp), 3), round(float(ss), 3)]
                             for r, pp, ss in zip(rows, p, s)}}
    print(json.dumps({k: (v if k not in ("ridge", "gp") else {kk: vv for kk, vv in v.items() if kk != "loo"})
                      for k, v in out.items()}, indent=1))
    json.dump(out, open(MODEL, "w"), indent=1)
    if "--apply" in sys.argv:
        sites = [json.loads(l) for l in open(os.path.join(D, "sites.jsonl")) if l.strip()]
        props = {json.loads(l)["name"]: json.loads(l) for l in open(os.path.join(D, "props.jsonl"))}
        sites = [s for s in sites if "error" not in s]
        res = {}
        for kind in ("ridge", "gp"):
            X = feats([r["eb_gfn2"] for r in rows], [r["q1"] for r in rows], [r["EA"] for r in rows], kind)
            mdl = (Ridge() if kind == "ridge" else GP()).fit(X, y)
            Xs = feats([s["E_b_gfn2"] for s in sites], [s["q_Ag1"] for s in sites],
                       [props[s["name"]]["EA"] for s in sites], kind)
            p, sg = mdl.predict(Xs)
            for s, pp, ss in zip(sites, p, sg):
                res.setdefault(s["name"], {})[kind] = [round(float(pp), 3), round(float(ss), 3)]
        json.dump(res, open(PRED, "w"), indent=1)
        print(f"predicted {len(res)}")


if __name__ == "__main__":
    main()
