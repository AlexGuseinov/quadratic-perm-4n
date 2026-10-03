"""Estimate the difficulty of a CNF by solving a few random cubes.

Usage: python3 cube_sample.py base.cnf "73,74,...,84" NCUBES TIMELIMIT SEED
Each cube fixes the listed variables to random values (unit clauses appended in a temporary
copy) and runs kissat with the time limit. Prints status and time per cube.
"""
import os, random, subprocess, sys, tempfile, time
HERE = os.path.dirname(os.path.abspath(__file__))
base, varlist, ncubes, tl, seed = sys.argv[1], [int(v) for v in sys.argv[2].split(",")], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
TMP = os.environ.get("SYM_TMP", tempfile.gettempdir())
with open(base) as f:
    lines = f.readlines()
hdr = next(i for i, l in enumerate(lines) if l.startswith("p cnf"))
nv, nc = map(int, lines[hdr].split()[2:4])
random.seed(seed)
for k in range(ncubes):
    units = [v if random.random() < 0.5 else -v for v in varlist]
    path = os.path.join(TMP, f"cube_{os.getpid()}_{k}.cnf")
    with open(path, "w") as g:
        g.writelines(lines[:hdr])
        g.write(f"p cnf {nv} {nc + len(units)}\n")
        g.writelines(lines[hdr + 1:])
        for u in units:
            g.write(f"{u} 0\n")
    t0 = time.time()
    r = subprocess.run([os.path.join(HERE, "..", "bin", "kissat"), "-q", f"--time={tl}", path], capture_output=True, text=True)
    st = {10: "SAT", 20: "UNSAT"}.get(r.returncode, "TIMEOUT")
    print(f"cube {k}: {st} {time.time() - t0:.1f}s", flush=True)
    os.remove(path)
