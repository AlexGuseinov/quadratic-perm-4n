"""Test Conjecture H on symmetric subspaces: for each symmetry case (sym_search.py)
with few free bits, count beta in the equivariant space whose 2^n - 1 components all
have rank n-2 (beta_enum mode 1). A nonzero count gives an n-dim constant-rank-(n-2)
subspace of Alt(n, F_2) (the components of such a beta span it, the map b -> B_b being
injective since no component vanishes).
Usage: python3 conjH_sym.py n MAXBITS [CASE ...]"""
import os, subprocess, sys, time, json
from itertools import combinations
from sym_search import cases, equivariant_basis
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")   # repository root (logs/ is there)
n, maxbits = int(sys.argv[1]), int(sys.argv[2])
cs = cases(n)
idxs = [int(x) for x in sys.argv[3:]] or range(len(cs))
npairs = n * (n - 1) // 2
for idx in idxs:
    p, la, A, lb, B = cs[idx]
    basis = equivariant_basis(A, B, n)[0]
    if len(basis) > maxbits:
        continue
    lines = [f"{n} {len(basis)} 1"] + ["".join(f"{(v >> (n * q)) & ((1 << n) - 1):02x}" for q in range(npairs)) for v in [0] + basis]
    t0 = time.time()
    r = subprocess.run([os.path.join(HERE, "beta_enum")], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True)
    out = r.stdout.split("\n")
    cnt = int(out[-2].split()[1])
    rec = {"n": n, "case": idx, "p": p, "A": la, "B": lb, "bits": len(basis), "const_rank_betas": cnt,
           "seconds": round(time.time() - t0, 1), "example": out[0] if cnt else None}
    with open(os.path.join(ROOT, "logs", f"conjH_sym_n{n}.jsonl"), "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"case {idx:3d} p={p:<3d} A={la:<12s} B={lb:<12s} bits={len(basis):3d} const-rank betas={cnt}", flush=True)
