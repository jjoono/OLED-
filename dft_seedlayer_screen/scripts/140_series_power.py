"""How many derivatives, how many films, and how precise must the closure
thickness be, for an interior optimum to be detectable -- and how often does
noise fake one?

    python scripts/140_series_power.py

The derivative series (package I) is designed around one prediction: if a good
seed binds hard and leaves silver metallic, the best member is an interior one,
not HATCN. That prediction is worth something only if the experiment can tell
an interior minimum from a monotone trend plus noise. This simulates the
experiment before it is run.

  truth H1   closure thickness d(x) = d0 + 4*Delta*(x - x*)^2, members evenly
             spaced on the electronic axis x in [0, 1], x* inside; Delta is how
             much thinner the optimum closes than the worse end
  truth H0   monotone, d(x) = d0 + S*(1 - x): the best member is at the end
  noise      each film's closure thickness ~ N(true, sigma); r films per member
  steps      ex-situ thickness series read the onset to the nearest step;
             in-situ resistance during deposition does not quantise

Decision rule: quadratic fit to the member means; an interior optimum is
declared when the curvature is positive at one-sided p < 0.05 and the vertex
lies strictly between the first and last member. Power is how often it is
declared under H1; the false-positive rate is how often under H0.
"""
import os, json, itertools

import numpy as np
from scipy import stats

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RNG = np.random.default_rng(7)
TRIALS = 4000


def measure(d, sigma, r, step):
    y = d[:, None] + sigma * RNG.standard_normal((len(d), r))
    if step:
        y = np.ceil(y / step) * step      # onset read at the next step up
    return y.mean(axis=1), y.std(axis=1, ddof=1) if r > 1 else None


def declare(x, ym, k_films):
    X = np.stack([np.ones_like(x), x, x ** 2], axis=1)
    beta, res, *_ = np.linalg.lstsq(X, ym, rcond=None)
    dof = len(x) - 3
    if dof < 1:
        return False
    s2 = float(np.sum((ym - X @ beta) ** 2)) / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    t = beta[2] / np.sqrt(max(cov[2, 2], 1e-30))
    p = 1 - stats.t.cdf(t, dof)
    if beta[2] <= 0 or p >= 0.05:
        return False
    v = -beta[1] / (2 * beta[2])
    return x[0] < v < x[-1]


def rate(k, sigma, r, step, delta=None, slope=None, xstar=0.55):
    x = np.linspace(0, 1, k)
    if delta is not None:
        # normalised so the worse end sits exactly delta above the optimum
        d = 7.0 + delta * (x - xstar) ** 2 / max(xstar, 1 - xstar) ** 2
    else:
        d = 7.0 + slope * (1 - x)
    hits = 0
    for _ in range(TRIALS):
        ym, _ = measure(d, sigma, r, step)
        hits += declare(x, ym, r)
    return hits / TRIALS


def main():
    out = {"power": [], "false_positive": []}
    print("POWER: P(interior optimum declared | it is real)")
    print(f"{'members':>8}{'films/mem':>10}{'sigma nm':>9}{'readout':>10}"
          f"{'Delta 0.3':>11}{'Delta 0.5':>11}{'Delta 1.0':>11}")
    for k, r, sigma, step in itertools.product((4, 5, 6), (1, 2, 3), (0.2, 0.3, 0.5),
                                               (None, 0.5)):
        pw = [rate(k, sigma, r, step, delta=D) for D in (0.3, 0.5, 1.0)]
        out["power"].append({"members": k, "films": r, "sigma": sigma,
                             "step": step, "power": dict(zip(("0.3", "0.5", "1.0"), pw))})
        print(f"{k:>8}{r:>10}{sigma:>9.1f}{('step 0.5' if step else 'in-situ'):>10}"
              + "".join(f"{p:>11.2f}" for p in pw))
    print("\nFALSE POSITIVES: P(interior optimum declared | truth is monotone, "
          "end member best by 1 nm)")
    for k, r, sigma, step in itertools.product((4, 5, 6), (1, 3), (0.3, 0.5), (None, 0.5)):
        fp = rate(k, sigma, r, step, slope=1.0)
        out["false_positive"].append({"members": k, "films": r, "sigma": sigma,
                                      "step": step, "rate": fp})
        print(f"{k:>8}{r:>10}{sigma:>9.1f}{('step 0.5' if step else 'in-situ'):>10}"
              f"{fp:>11.3f}")
    json.dump(out, open(os.path.join(ROOT, "runs", "series_power.json"), "w"), indent=1)
    print("\nwrote runs/series_power.json")


if __name__ == "__main__":
    main()
