r"""Run the campaign: several folders at a time, each folder in path order.

    python run_campaign.py [--workers N] [--g16 C:\G16W]

Each candidate folder is one orbital chain with its own .chk, so folders are
the unit of parallelism and no chain has to be split. One worker owns a folder
from start to finish and runs its jobs in the order ORDER.txt gives, which is
not the alphabetical order: the reference re-run <name>_t0_z0_ref must go last,
after the chain it is meant to read.

Each worker gets its own Gaussian scratch directory. Sharing one is what stops
several G16W processes from running at the same time.

A job whose .out already ends in "Normal termination" is skipped, so the script
can be stopped and restarted and will pick up where it left off.
"""
import argparse, os, shutil, subprocess, sys, threading, time

HERE = os.path.dirname(os.path.abspath(__file__))
LOCK = threading.Lock()
STATE = {"done": 0, "skip": 0, "bad": 0}


def finished(out):
    """True if this job already ran to completion."""
    try:
        with open(out, "rb") as f:
            f.seek(max(0, os.path.getsize(out) - 4096))
            return b"Normal termination" in f.read()
    except OSError:
        return False


def say(msg):
    with LOCK:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def jobs_of(folder):
    """Job stems in the order they must run."""
    stems = sorted(f[:-4] for f in os.listdir(folder) if f.endswith(".gjf"))
    order = os.path.join(folder, "ORDER.txt")
    if os.path.exists(order):
        names = [l.strip() for l in open(order) if l.strip()]
        listed = [n for n in names if n + ".gjf" in os.listdir(folder)]
        # Anything present but not listed is a job added after the folder was
        # written -- a top-up height, say. Run it after what the file names
        # rather than ignoring it, which is what filtering to the list alone
        # would do.
        extra = [n for n in stems if n not in set(listed)]
        return listed + [n for n in extra if not n.endswith("_ref")] \
            + [n for n in extra if n.endswith("_ref")]
    # No ORDER.txt: fall back to alphabetical but push the reference re-runs to
    # the end, which is the one thing alphabetical order gets wrong.
    return ([s for s in stems if not s.endswith("_ref")]
            + [s for s in stems if s.endswith("_ref")])


def _threads_of(folder):
    """%NProcShared this folder's jobs ask for."""
    for fn in sorted(os.listdir(folder)):
        if fn.endswith(".gjf"):
            for line in open(os.path.join(folder, fn)):
                if line.lower().startswith("%nprocshared="):
                    try:
                        return int(line.split("=")[1])
                    except ValueError:
                        break
            break
    return 8


def _cost(folder):
    """Rough work in a folder: atoms cubed times the number of jobs."""
    gjfs = [f for f in os.listdir(folder) if f.endswith(".gjf")]
    if not gjfs:
        return 0.0
    atoms = sum(1 for line in open(os.path.join(folder, gjfs[0]))
                if len(line.split()) == 4 and line.split()[0].isalpha())
    return (atoms ** 3) * len(gjfs)


def run_folder(folder, g16, scr):
    os.makedirs(scr, exist_ok=True)
    env = dict(os.environ)
    env["GAUSS_EXEDIR"] = g16
    env["GAUSS_SCRDIR"] = scr
    env["PATH"] = g16 + os.pathsep + env.get("PATH", "")
    name = os.path.basename(folder)
    exe = os.path.join(g16, "g16.exe")

    for stem in jobs_of(folder):
        out = os.path.join(folder, stem + ".out")
        if finished(out):
            with LOCK:
                STATE["skip"] += 1
            continue
        say(f"{name}: {stem}")
        t0 = time.time()
        try:
            subprocess.run([exe, stem + ".gjf", stem + ".out"],
                           cwd=folder, env=env, check=False)
        except OSError as e:
            say(f"{name}: cannot start Gaussian -- {e}")
            with LOCK:
                STATE["bad"] += 1
            return
        dt = time.time() - t0
        ok = finished(out)
        with LOCK:
            STATE["done" if ok else "bad"] += 1
        if not ok:
            # Do not go on down a chain whose orbitals were never written: every
            # later point would read a stale .chk and the barrier would be made
            # of two different electronic states.
            say(f"{name}: {stem} did NOT finish normally after {dt/60:.1f} min "
                f"-- stopping this folder")
            return
        say(f"{name}: {stem} done in {dt/60:.1f} min")
        for junk in os.listdir(folder):
            if junk.startswith("Gau-") or junk == "fort.7":
                try:
                    os.remove(os.path.join(folder, junk))
                except OSError:
                    pass
    say(f"{name}: folder complete")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=0,
                    help="folders to run at once (default: cores // 8)")
    ap.add_argument("--g16", default=r"C:\G16W")
    ap.add_argument("--threads", type=int, default=0,
                    help="physical cores to keep busy (default: cores // 2)")
    ap.add_argument("--scratch", default="")
    a = ap.parse_args()

    if not os.path.exists(os.path.join(a.g16, "g16.exe")):
        print(f"Gaussian not found at {os.path.join(a.g16, 'g16.exe')}")
        print("Pass the right one:  python run_campaign.py --g16 C:\\path\\to\\G16W")
        return 1

    cores = os.cpu_count() or 8
    # os.cpu_count() reports logical processors, and SMT threads do not help a
    # Fock build -- they contend for the same floating-point units. Budget on
    # physical cores, assuming SMT is on, and give each job the 8 threads its
    # %NProcShared asks for.
    workers = a.workers or max(1, (cores // 2) // 8)
    budget = a.threads or max(8, cores // 2)
    scratch = a.scratch or os.path.join(
        os.environ.get("TEMP", HERE), "gauscr")

    folders = [os.path.join(HERE, d) for d in os.listdir(HERE)
               if os.path.isdir(os.path.join(HERE, d))
               and any(f.endswith(".gjf") for f in os.listdir(os.path.join(HERE, d)))]
    # Longest first. Folders taken alphabetically leave the 85-atom Bphen dimer
    # grinding on one worker for a day after everything else has finished;
    # starting the long poles immediately overlaps them with the short ones,
    # which is the whole of the difference in when the campaign ends.
    folders.sort(key=lambda f: -_cost(f))
    if not folders:
        print("no candidate folders next to this script")
        return 1

    print(f"{len(folders)} folders, {budget} physical cores to fill, "
          f"{cores} logical processors visible")
    print("  order: " + ", ".join(os.path.basename(f) for f in folders[:5])
          + ", ...")
    print(f"scratch under {scratch}\n")

    # Folders are claimed against a thread budget, not a worker count, because
    # they no longer all ask for the same number of threads: an 85-atom dimer
    # gets 16 where a 20-atom complex gets 8. Counting folders would either
    # oversubscribe the machine or leave half of it idle.
    queue = list(folders)
    qlock = threading.Condition()
    free = [budget]

    def worker(k):
        scr = os.path.join(scratch, f"w{k}")
        while True:
            with qlock:
                while True:
                    if not queue:
                        return
                    i = next((j for j, f in enumerate(queue)
                              if _threads_of(f) <= free[0]), None)
                    if i is not None:
                        break
                    qlock.wait(30)
                folder = queue.pop(i)
                want = _threads_of(folder)
                free[0] -= want
            try:
                run_folder(folder, a.g16, scr)
            finally:
                with qlock:
                    free[0] += want
                    qlock.notify_all()

    t0 = time.time()
    nthread = max(workers, budget // 8)
    threads = [threading.Thread(target=worker, args=(k,), daemon=True)
               for k in range(nthread)]
    for t in threads:
        t.start()
    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        print("\nstopped -- re-run this script to continue where it left off")
        return 1

    print(f"\n{STATE['done']} ran, {STATE['skip']} already done, "
          f"{STATE['bad']} failed, {(time.time()-t0)/3600:.1f} h wall")
    if STATE["bad"]:
        print("Re-run this script to retry the folders that stopped.")
    shutil.rmtree(scratch, ignore_errors=True)
    print("When everything is done, send the whole gaussian_jobs folder back.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
