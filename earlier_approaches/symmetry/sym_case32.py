"""Case p = 3, A = B = f0+f0+f0 (n = 8), split by a conjugacy class in GL(3, 4).

After the normal form (sym_search.normal_form): V1 = <e6, e7> is a kernel line and
phi_x = beta(x, .) restricted to V_w = <e0..e5> is, for x in V1 \\ 0, an F_4-linear
bijection V_w -> W_w (A = B = omega on these parts; N5 and m_A = m_B), with phi_{e6}
the identity (SB2''). Pairs (g, g) with g in GL_{F_4}(V_w), acting on input and output,
keep all this and replace M = phi_{e6}^{-1} phi_{e7} by g M g^{-1}. M and M + I are
invertible (phi_{e7}, phi_{e6+e7} are bijective). So WLOG M is one representative per
conjugacy class of such matrices; then beta(e7, e_{2i}) = M e_{2i} is fixed (i = 0, 1, 2)
and each class is decided by kissat on the full model (with a DRAT check when short).
"""
import json, os, time
from itertools import combinations, product
import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from sym_search import cases, normal_form, affine_space
import subprocess, tempfile
from sym_search import build_sym
from sym_runner import write_cnf, KISSAT, DRAT
from du4perm_sat import extract, lut, check

# F_4 = {0, 1, w, w^2} as 2-bit ints: bit0 = coefficient of 1, bit1 = coefficient of w
def f4mul(a, b):
    r = 0
    # (a0 + a1 w)(b0 + b1 w) = a0b0 + (a0b1 + a1b0) w + a1b1 w^2, w^2 = w + 1
    a0, a1, b0, b1 = a & 1, a >> 1, b & 1, b >> 1
    c0 = (a0 & b0) ^ (a1 & b1)
    c1 = (a0 & b1) ^ (a1 & b0) ^ (a1 & b1)
    return c0 | (c1 << 1)

def mmul(X, Y):
    return tuple(tuple(
        (lambda s: s)(0 if False else __import__("functools").reduce(lambda acc, k: acc ^ f4mul(X[i][k], Y[k][j]), range(3), 0))
        for j in range(3)) for i in range(3))

def det(X):
    m = f4mul
    return (m(X[0][0], m(X[1][1], X[2][2]) ^ m(X[1][2], X[2][1]))
            ^ m(X[0][1], m(X[1][0], X[2][2]) ^ m(X[1][2], X[2][0]))
            ^ m(X[0][2], m(X[1][0], X[2][1]) ^ m(X[1][1], X[2][0])))

INV = {1: 1, 2: 3, 3: 2}
def inv(X):
    d = INV[det(X)]
    m = f4mul
    cof = [[0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            r = [x for x in range(3) if x != i]; c = [y for y in range(3) if y != j]
            minor = m(X[r[0]][c[0]], X[r[1]][c[1]]) ^ m(X[r[0]][c[1]], X[r[1]][c[0]])
            cof[j][i] = m(d, minor)          # adjugate (signs are +1 in char 2)
    return tuple(tuple(row) for row in cof)

I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
allM = [tuple(tuple(v[3 * i:3 * i + 3]) for i in range(3)) for v in product(range(4), repeat=9)]
good = [M for M in allM if det(M) and det(tuple(tuple(M[i][j] ^ I3[i][j] for j in range(3)) for i in range(3)))]
# GL(3,4) generators: a diagonal matrix with w and an elementary + permutation matrix
G1 = ((2, 0, 0), (0, 1, 0), (0, 0, 1))
G2 = ((1, 1, 0), (0, 1, 0), (0, 0, 1))
G3 = ((0, 0, 1), (1, 0, 0), (0, 1, 0))
gens = [(g, inv(g)) for g in (G1, G2, G3)]
seen, reps = set(), []
goodset = set(good)
for M in good:
    if M in seen:
        continue
    reps.append(M)
    stack = [M]; seen.add(M)
    while stack:
        X = stack.pop()
        for g, gi in gens:
            Y = mmul(mmul(g, X), gi)
            if Y not in seen:
                seen.add(Y); stack.append(Y)
assert seen == goodset
print(f"{len(good)} matrices, {len(reps)} conjugacy classes", flush=True)

n = 8
p, la, A, lb, B = cases(n)[32]
ok, why, base = normal_form(p, A, B, n, la, lb)
pairs = list(combinations(range(n), 2)); pidx = {pq: t for t, pq in enumerate(pairs)}
def f4vec_to_bits(col):     # F_4^3 column -> bits on e0..e5 (e_{2i} = 1, e_{2i+1} = w)
    v = 0
    for i, c in enumerate(col):
        v |= (c & 1) << (2 * i) | (c >> 1) << (2 * i + 1)
    return v
log = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "logs", "sym_case32.jsonl")
TMP = os.environ.get("SYM_TMP", tempfile.gettempdir())
import sys
TL = int(sys.argv[1]) if len(sys.argv) > 1 else 3600
for k, M in enumerate(reps):
    rows = list(base)
    for i in range(3):
        target = f4vec_to_bits([M[r][i] for r in range(3)])   # M e_{2i}
        for t in range(n):
            rows.append((1 << (pidx[(2 * i, 7)] * n + t), target >> t & 1))
    t0 = time.time()
    rec = {"class": k, "M": M}
    try:
        enc, c, ell, dim = build_sym(n, A, B, True, False, rows)
    except ValueError as e:
        rec.update(status="UNSAT", note=str(e))
    else:
        cnf = os.path.join(TMP, f"case32_{k}_{os.getpid()}.cnf")
        write_cnf(enc, cnf, f"case 32 class {k}")
        r = subprocess.run([KISSAT, "-q", f"--time={TL}", cnf], capture_output=True, text=True)
        rec.update(bits=dim, vars=enc.pool.top, clauses=len(enc.clauses))
        if r.returncode == 10:
            model = [int(x) for line in r.stdout.splitlines() if line.startswith("v ")
                     for x in line[2:].split() if x != "0"]
            C, L = extract(model, c, ell, n)
            rec.update(status="SAT", check=check(lut(C, L, n), n), lut=lut(C, L, n))
        elif r.returncode == 20:
            rec["status"] = "UNSAT"
            if time.time() - t0 < 300:
                proof = cnf + ".drat"
                subprocess.run([KISSAT, "-q", "--no-binary", cnf, proof], capture_output=True)
                v = subprocess.run([DRAT, cnf, proof], capture_output=True, text=True)
                rec["drat"] = "VERIFIED" if "s VERIFIED" in v.stdout.replace("\r", "\n") else "NOT VERIFIED"
                os.remove(proof)
        else:
            rec["status"] = "TIMEOUT"
        os.remove(cnf)
    rec.update(seconds=round(time.time() - t0, 1), when=time.strftime("%Y-%m-%d %H:%M:%S"))
    with open(log, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print({k2: v2 for k2, v2 in rec.items() if k2 != "lut"}, flush=True)
