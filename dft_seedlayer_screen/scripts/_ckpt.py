"""Checkpoint helper: make a calculation survive the container being replaced.

THE PROBLEM. This container does not merely kill processes -- it is restored from
an earlier snapshot. Evidence accumulated over six occurrences: the local .git
reverts to an OLD commit (a fresh clone would give the remote HEAD), the working
tree reverts with it, pip installs vanish, /tmp is wiped, and setsid/nohup/disown
make no difference. Nothing running inside can survive, because the machine state
is replaced rather than signalled.

So there is no process-level fix. The only thing that survives is the git remote,
and the only useful question is how much work a rollback costs. This module makes
that answer "one unit", by writing each completed unit to a JSON and pushing it
immediately.

USAGE

    from _ckpt import Checkpoint
    ck = Checkpoint("runs/myresults.json")

    for job in jobs:
        if ck.has(job):            # already done in an earlier container
            continue
        value = expensive(job)
        ck.put(job, value)         # writes, commits, pushes

WHY THE SCRIPT AND NOT A WATCHER PROCESS. A separate guard loop works until the
rollback takes the guard too -- and it always does, since it lives in the same
container. Putting the checkpoint inside the calculation means the push happens
in the same breath as the result.

GRANULARITY. Make the unit small enough that redoing one is cheap. A six-point
scan checkpointed per COMPLEX loses the whole scan; checkpointed per SCAN POINT
it loses one SCF. The cost of a push is ~1 s against SCF minutes, so prefer the
finer unit.

ON IDLE. The rollbacks correlate with the session going idle while a background
job runs. Prefer running work in FOREGROUND chunks inside a turn over backgrounding
it and waiting -- the harness caps a call at ten minutes, so size units to fit.
"""
import json, os, subprocess, time

REPO = "/home/user/OLED-"
BRANCH = "claude/project-setup-conda-ftxhnl"


def _git(*args, check=False):
    return subprocess.run(["git", "-C", REPO, *args], capture_output=True,
                          text=True, check=check)


class Checkpoint:
    def __init__(self, path, push=True, label="checkpoint"):
        """path: repo-relative or absolute path to the results JSON."""
        self.abspath = path if os.path.isabs(path) else os.path.join(REPO, path)
        self.relpath = os.path.relpath(self.abspath, REPO)
        self.push = push
        self.label = label
        os.makedirs(os.path.dirname(self.abspath), exist_ok=True)
        self.data = self._load()

    def _load(self):
        try:
            with open(self.abspath) as f:
                return json.load(f)
        except Exception:
            return {}

    def has(self, key):
        """True if this unit already has a real result (errors do not count)."""
        v = self.data.get(key)
        return v is not None and not (isinstance(v, dict) and "error" in v)

    def get(self, key, default=None):
        return self.data.get(key, default)

    def put(self, key, value, commit=True):
        self.data[key] = value
        # Re-read first: a concurrent writer (or a resumed run) may hold units this
        # process never computed, and clobbering them would silently lose work.
        merged = self._load()
        merged.update(self.data)
        self.data = merged
        tmp = self.abspath + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.data, f, indent=2, sort_keys=True)
        os.replace(tmp, self.abspath)      # atomic: a rollback mid-write cannot
                                           # leave a truncated JSON behind
        if commit:
            self.commit(f"{self.label}: {key}")

    def commit(self, message):
        """Commit and push. Never raises -- a push failure must not lose the run."""
        try:
            _git("add", self.relpath)
            st = _git("diff", "--cached", "--quiet")
            if st.returncode == 0:
                return False                      # nothing staged
            body = (f"{message}\n\n"
                    "Auto-checkpoint. The container is restored from snapshots, so\n"
                    "only what reaches the remote survives; this pins one completed\n"
                    "unit of a long calculation.\n\n"
                    "Co-Authored-By: Claude Opus 5 <noreply@anthropic.com>\n"
                    "Claude-Session: https://claude.ai/code/session_01QdGePAvo7ESgQUdueHRqPN")
            _git("commit", "-q", "-m", body)
            if self.push:
                for attempt in range(4):
                    r = _git("push", "-q", "origin", BRANCH)
                    if r.returncode == 0:
                        return True
                    # diverged: another container pushed. Rebase and retry rather
                    # than force -- the other side's units are results too.
                    _git("fetch", "-q", "origin", BRANCH)
                    _git("rebase", "-q", f"origin/{BRANCH}")
                    time.sleep(2 ** attempt)
                print("  [ckpt] push failed after retries; result is committed "
                      "locally only and will not survive a rollback", flush=True)
            return True
        except Exception as exc:
            print(f"  [ckpt] commit failed ({type(exc).__name__}); continuing",
                  flush=True)
            return False


def resume_banner(ck, jobs):
    done = [j for j in jobs if ck.has(j)]
    todo = [j for j in jobs if not ck.has(j)]
    print(f"[ckpt] {len(done)}/{len(jobs)} units already done"
          f"{' -> ' + ', '.join(done) if done else ''}", flush=True)
    if todo:
        print(f"[ckpt] remaining: {', '.join(todo)}", flush=True)
    return todo
