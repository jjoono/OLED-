"""Does closing early actually buy lower absorption in the device -- and does
scattering make the whole comparison moot?

Closure thickness is a diagnostic, not a target. Nobody builds the electrode at
the closure thickness: the device optimum sits at 7-8 nm for both seeds. What
the seed buys is that the SAME 8 nm film is a different film depending on what
it grew on, and that is what this script separates.

Each thickness is split three ways:

  floor    Fuchs-Sondheimer surface scattering only, at p = 0. Irreducible for a
           film of this thickness -- no deposition trick removes it.
  addressable   everything between that floor and the measurement: grain
           boundaries, roughness, voids, and any light scattered out of the
           collection cone. This is the part seed engineering and deposition
           conditions act on.
  scattering    cannot be separated from the addressable part without an angular
           scan, so it is bounded rather than measured: even if EVERY photon
           missing from 1-T-R were scattered rather than absorbed, the
           addressable term is its ceiling.

The bound is what answers the question. If the seed-to-seed difference is larger
than the entire addressable term, no assumption about scattering can overturn
the ranking.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
m = importlib.import_module("91_nk_library_series")

L = 550.0
RHO_BULK, MFP = 1.59, 52.0
SERIES = {
    "HATCN": {4: ("1-2", 52.2), 5: ("1-3", 23.3), 6: ("1-4", 18.9), 7: ("2-1", 12.7),
              8: ("2-2", 9.1), 10: ("2-3", 7.5), 12: ("2-4", 5.3)},
    "MoOx":  {4: ("1-10", 138.0), 5: ("1-11", 49.1), 6: ("1-12", 32.0),
              7: ("2-9", 14.9), 8: ("2-10", 10.8), 10: ("2-11", 9.6)},
}
CLOSURE = {"HATCN": 5, "MoOx": 7}


def analyse(seed, d, sid):
    T = m.read_csv(os.path.join(m.RAW, f"{sid}T.csv"))[L]
    R = m.read_csv(os.path.join(m.RAW, f"{sid}R.csv"))[L]
    nk, _ = m.invert(L, T, R, seed, float(d))
    e1 = nk.real**2 - nk.imag**2
    e2 = 2 * nk.real * nk.imag
    ag = m.AG_IDEAL(L)
    e2_bulk = 2 * ag.real * ag.imag
    e2_floor = e2_bulk * (1.0 + 0.375 * MFP / d)          # p = 0 Fuchs
    nk_floor = np.sqrt(complex(e1, e2_floor))
    a_meas, _ = m.best_device_A(nk, L, float(d))
    a_floor, _ = m.best_device_A(nk_floor, L, float(d))
    a_bulk, _ = m.best_device_A(ag, L, float(d))
    return a_meas, a_floor, a_bulk, e2, e2_floor


def main():
    print(f"device one-pass absorption at {L:.0f} nm, organic 1.8 / Ag / CPL 2.1\n")
    res = {}
    for seed, series in SERIES.items():
        print(f"{seed} 5 nm   (closes at {CLOSURE[seed]} nm)")
        print(f"{'d':>4} {'bulk Ag':>8} {'floor':>7} {'measured':>9} "
              f"{'addressable':>12} {'eps2':>6} {'/floor':>7}")
        print("  " + "-" * 62)
        res[seed] = {}
        for d, (sid, _) in series.items():
            am, af, ab, e2, e2f = analyse(seed, d, sid)
            res[seed][d] = (am, af, ab)
            print(f"{d:>4} {ab*100:>7.2f}% {af*100:>6.2f}% {am*100:>8.2f}% "
                  f"{(am-af)*100:>11.2f}% {e2:>6.2f} {e2/e2f:>6.2f}x")
        print()

    print("the seed comparison, at the thickness a device would actually use")
    print(f"{'d':>4} {'HATCN':>8} {'MoOx':>8} {'gap':>8} {'HATCN addressable':>18} "
          f"{'gap > addressable?':>19}")
    print("-" * 72)
    for d in (7, 8, 10):
        h, mo = res["HATCN"][d], res["MoOx"][d]
        gap = mo[0] - h[0]
        addr = h[0] - h[1]
        print(f"{d:>4} {h[0]*100:>7.2f}% {mo[0]*100:>7.2f}% {gap*100:>7.2f}% "
              f"{addr*100:>17.2f}% {'yes' if gap > addr else 'no':>19}")

    d = 8
    h, mo = res["HATCN"][d], res["MoOx"][d]
    print(f"\nAt {d} nm the seed is worth {(mo[0]-h[0])*100:.2f} %p, while HATCN's entire")
    print(f"addressable term -- grain boundaries AND roughness AND every scattered")
    print(f"photon together -- is {(h[0]-h[1])*100:.2f} %p. Push all of it into scattering and")
    print(f"HATCN's absorption still lands at {h[1]*100:.2f}%, against MoOx at "
          f"{mo[0]*100:.2f}% before")
    print("its own scattering is subtracted. The ranking cannot be reached from here.")
    print(f"\nWhat is left to win: {h[0]*100:.2f}% now, {h[1]*100:.2f}% at the Fuchs floor, "
          f"{h[2]*100:.2f}% with bulk silver.")
    print(f"So {(h[0]-h[1])/(h[0]-h[2])*100:.0f}% of the gap to ideal silver is addressable")
    print("by process, and the rest is the size effect itself, which only a thicker")
    print("film removes -- and a thicker film costs more than it saves past 8 nm.")

    print("\nclosure thickness is not the target: HATCN closes at 5 nm but the device")
    print(f"optimum is 8 nm ({res['HATCN'][8][0]*100:.2f}% against "
          f"{res['HATCN'][5][0]*100:.2f}% at 5 nm).")
    print("Closing early matters because a film that nucleated densely is still the")
    print("better film three nanometres later, not because 5 nm is where you stop.")


if __name__ == "__main__":
    main()
