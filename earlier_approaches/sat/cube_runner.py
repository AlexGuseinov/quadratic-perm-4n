"""Checkpointed cube-and-conquer runner for the n = 8 constant-rank SAT problems.

The search space is split into 2^k cubes by fixing k basis variables of B2
(B1 is already fixed to the canonical rank-6 form). Each cube is solved with a
conflict budget; its result (UNSAT / SAT / TIMEOUT) is appended to a JSON-lines
progress file, so work survives machine restarts. Timed-out cubes are later
split further (depth increases) and retried.

Usage:
  python3 cube_runner.py MODE WORKER NWORKERS [--k K] [--conflicts C]
  MODE: du4 (constant rank 6 + s_a >= 2) or H (constant rank 6 only)
Progress: logs/cubes_<MODE>.jsonl ; a SAT solution is written to
logs/SOLUTION_<MODE>.txt and all workers stop.
"""

import json
import os
import sys
import time

from pysat.solvers import Solver

import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from constrank_sat import build, decode, verify, s_profile

N, R, D = 8, 6, 8
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")   # repository root (logs/ is there)


def load_progress(path):
    done = {}
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    rec = json.loads(line)
                    done[rec["cube"]] = rec["status"]
    return done


def cube_key(bits):
    return "".join("1" if b else "0" for b in bits)


def children(key, extra):
    out = [key]
    for _ in range(extra):
        out = [k + "0" for k in out] + [k + "1" for k in out]
    return out


def pending_cubes(done, k, split_extra):
    """Top-level cubes of length k; timed-out cubes are refined by split_extra
    more variables. A cube is resolved if it or all its refinements are UNSAT."""
    todo = []

    def visit(key):
        st = done.get(key)
        if st == "UNSAT":
            return
        if st == "TIMEOUT":
            for ch in children(key, split_extra)[0:]:
                if ch != key:
                    visit(ch)
            return
        if st is None:
            todo.append(key)

    for i in range(1 << k):
        visit(format(i, f"0{k}b"))
    return todo


def main():
    mode = sys.argv[1]
    worker, nworkers = int(sys.argv[2]), int(sys.argv[3])
    k = int(sys.argv[sys.argv.index("--k") + 1]) if "--k" in sys.argv else 10
    # conflict budget per cube (CaDiCaL ignores pysat's interrupt; ~1.7k conflicts/s here)
    budget = int(sys.argv[sys.argv.index("--conflicts") + 1]) if "--conflicts" in sys.argv else 500000
    split_extra = 3
    du4 = mode == "du4"
    prog = os.path.join(ROOT, "logs", f"cubes_{mode}.jsonl")
    sol_path = os.path.join(ROOT, "logs", f"SOLUTION_{mode}.txt")

    enc, basis, pairs = build(N, R, D, du4=du4)
    split_vars = basis[1] + basis[2]  # B2 then B3 entries, 56 variables
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        while True:
            if os.path.exists(sol_path):
                print("solution already found; stopping", flush=True)
                return
            done = load_progress(prog)
            # a worker owns a top-level cube and all its refinements
            todo = [c for c in pending_cubes(done, k, split_extra)
                    if int(c[:k], 2) % nworkers == worker]
            if not todo:
                print("no pending cubes for this worker", flush=True)
                return
            key = todo[0]
            assumptions = [split_vars[i] if ch == "1" else -split_vars[i]
                           for i, ch in enumerate(key)]
            t0 = time.time()
            s.conf_budget(budget)
            res = s.solve_limited(assumptions=assumptions)
            dt = time.time() - t0
            status = {True: "SAT", False: "UNSAT", None: "TIMEOUT"}[res]
            rec = {"cube": key, "status": status, "seconds": round(dt, 1),
                   "worker": worker, "time": time.strftime("%Y-%m-%d %H:%M:%S")}
            with open(prog, "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(rec, flush=True)
            if res is True:
                mats = decode(s.get_model(), basis, pairs, N)
                with open(sol_path, "w") as f:
                    f.write(f"mode={mode} cube={key}\n")
                    f.write(f"rank check: {verify(mats, N, R)}\n")
                    f.write(f"s_a histogram: {s_profile(mats, N)}\n")
                    for i, rows in enumerate(mats):
                        f.write(f"B{i + 1}: " + " ".join(
                            f"{x:0{N}b}"[::-1] for x in rows) + "\n")
                print("SAT — solution saved", flush=True)
                return


if __name__ == "__main__":
    main()
