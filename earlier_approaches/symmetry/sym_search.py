"""Quadratic DU-4 permutations of F_2^n whose quadratic part has a prescribed
linear symmetry:   beta(A x, A y) = B beta(x, y)   for all x, y,
with A of prime order p and B^p = I (B = I allowed).

Why this covers every nontrivial symmetry.  If F(Ax + c) = B F(x) + (affine)
with (A, B) != (I, I), comparing quadratic parts gives beta(Ax, Ay) = B beta(x, y).
A = I forces B = I (beta spans F_2^n), so A != I; replacing (A, B) by a power we
may assume A has prime order p, and then B^p = I because beta spans F_2^n.
Up to conjugation (input basis) A is one representative per conjugacy class of
cyclic subgroups of order p, and B (output basis) runs over all conjugacy
classes of elements with B^p = I.

The linear part F(e_i) is left free (no symmetry assumed on it).
Usage:
  python3 sym_search.py list n                  # list the cases
  python3 sym_search.py cnf n CASE out.cnf      # write DIMACS for one case
  python3 sym_search.py solve n CASE            # solve with CaDiCaL (small n)
"""
import sys
import time
from itertools import combinations, combinations_with_replacement

import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
from du4perm_sat import Enc, T, F_, add_core, extract, lut, check


# ---------- GF(2)[x] helpers (polynomials as ints, bit i = coeff of x^i)
def deg(a):
    return a.bit_length() - 1


def pmod(a, m):
    dm = deg(m)
    while a and deg(a) >= dm:
        a ^= m << (deg(a) - dm)
    return a


def pmulmod(a, b, m):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> deg(m) & 1:
            a ^= m
    return r


def ppowmod(a, e, m):
    r = 1
    a = pmod(a, m)
    while e:
        if e & 1:
            r = pmulmod(r, a, m)
        a = pmulmod(a, a, m)
        e >>= 1
    return r


def irreducible(f):
    d = deg(f)
    for g in range(2, 1 << (d // 2 + 1)):
        if deg(g) >= 1 and pmod(f, g) == 0:
            return False
    return True


def ord2(p):
    k, v = 1, 2 % p
    while v != 1:
        v = v * 2 % p
        k += 1
    return k


def cyclotomic_factors(p):
    """Irreducible factors of (x^p - 1)/(x - 1) over F_2 (p odd prime)."""
    d = ord2(p)
    return [f for f in range(1 << d, 1 << (d + 1))
            if irreducible(f) and ppowmod(2, p, f) == 1]


def root_power_map(factors, p):
    """perm[j][i] = index of the minimal polynomial of alpha^j, alpha a root of factors[i]."""
    perm = {}
    for j in range(1, p):
        row = []
        for f in factors:
            r = ppowmod(2, j, f)          # alpha^j as a polynomial in alpha
            for gi, g in enumerate(factors):
                # evaluate g(r) mod f
                acc, pw = 0, 1
                for k in range(deg(g) + 1):
                    if g >> k & 1:
                        acc ^= pw
                    pw = pmulmod(pw, r, f)
                if acc == 0:
                    row.append(gi)
                    break
        perm[j] = row
    return perm


# ---------- matrices as lists of column bitmasks
def companion(f):
    d = deg(f)
    return [1 << (i + 1) for i in range(d - 1)] + [f ^ (1 << d)]


def blockdiag(blocks, n):
    cols, off = [], 0
    for blk in blocks:
        for col in blk:
            cols.append(col << off)
        off += len(blk)
    while len(cols) < n:
        cols.append(1 << len(cols))
    assert len(cols) == n
    return cols


def apply(M, v):
    r = 0
    i = 0
    while v:
        if v & 1:
            r ^= M[i]
        v >>= 1
        i += 1
    return r


def matmul(M, N):
    return [apply(M, col) for col in N]


def identity(n):
    return [1 << i for i in range(n)]


def order(M, n, cap=300):
    P = M[:]
    for k in range(1, cap):
        if P == identity(n):
            return k
        P = matmul(M, P)
    return None


# ---------- conjugacy classes
def element_types(p, n):
    """All conjugacy classes of elements of order p in GL(n, 2): (label, matrix)."""
    out = []
    if p == 2:
        for k in range(1, n // 2 + 1):
            out.append((f"J2^{k}", blockdiag([[1, 3]] * k, n)))
        return out
    fac = cyclotomic_factors(p)
    d = deg(fac[0])
    for m in range(1, n // d + 1):
        for combo in combinations_with_replacement(range(len(fac)), m):
            lab = "+".join(f"f{i}" for i in combo)
            out.append((lab, blockdiag([companion(fac[i]) for i in combo], n)))
    return out


def group_types(p, n):
    """One generator per conjugacy class of cyclic subgroups of order p."""
    if p == 2:
        return element_types(2, n)
    fac = cyclotomic_factors(p)
    d = deg(fac[0])
    perm = root_power_map(fac, p)
    seen, out = set(), []
    for m in range(1, n // d + 1):
        for combo in combinations_with_replacement(range(len(fac)), m):
            orbit = {tuple(sorted(perm[j][i] for i in combo)) for j in range(1, p)}
            key = min(orbit)
            if key in seen:
                continue
            seen.add(key)
            lab = "+".join(f"f{i}" for i in key)
            out.append((lab, blockdiag([companion(fac[i]) for i in key], n)))
    return out


def primes_for(n):
    size = 1
    for i in range(n):
        size *= (1 << n) - (1 << i)
    return [q for q in range(2, 1 << n) if size % q == 0 and all(q % r for r in range(2, q))]


def cases(n):
    out = []
    for p in primes_for(n):
        for la, A in group_types(p, n):
            assert order(A, n) == p
            for lb, B in element_types(p, n) + [("I", identity(n))]:
                out.append((p, la, A, lb, B))
    return out


# ---------- equivariant beta: nullspace of the linear conditions
def equivariant_basis(A, B, n):
    pairs = list(combinations(range(n), 2))
    pidx = {pq: t for t, pq in enumerate(pairs)}
    N = len(pairs) * n

    def u(pq, t):
        return pidx[pq] * n + t

    rows = []
    for (i, j) in pairs:
        for t in range(n):
            row = 0
            for (k, l) in pairs:
                coef = ((A[i] >> k & 1) & (A[j] >> l & 1)) ^ ((A[i] >> l & 1) & (A[j] >> k & 1))
                if coef:
                    row ^= 1 << u((k, l), t)
            for s in range(n):
                if B[s] >> t & 1:
                    row ^= 1 << u((i, j), s)
            if row:
                rows.append(row)
    # reduced row echelon form
    piv = {}
    for r in rows:
        for col, pr in piv.items():
            if r >> col & 1:
                r ^= pr
        if r:
            col = r.bit_length() - 1
            for c2 in list(piv):
                if piv[c2] >> col & 1:
                    piv[c2] ^= r
            piv[col] = r
    free = [col for col in range(N) if col not in piv]
    basis = []
    for fcol in free:
        v = 1 << fcol
        for col, pr in piv.items():
            if pr >> fcol & 1:
                v |= 1 << col
        basis.append(v)
    return basis, pairs, u


def trailing_fixed(M, n):
    """Number of trailing basis vectors fixed by M (the identity tail of blockdiag)."""
    f = 0
    for i in range(n - 1, -1, -1):
        if M[i] != 1 << i:
            break
        f += 1
    return f


def fixed_dim(M, n):
    """dim ker(M - I)."""
    rows = [M[i] ^ (1 << i) for i in range(n)]
    rank = 0
    for bit in range(n):
        piv = next((r for r in rows if r >> bit & 1), None)
        if piv is None:
            continue
        rows = [r ^ piv if (r >> bit & 1) and r is not piv else r for r in rows if r is not piv]
        rank += 1
    return n - rank


def normal_form(p, A, B, n, la=None, lb=None):
    """Extra WLOG conditions for odd p (proofs in NOTES_lemmas.md, "Symmetric case").

    Returns (ok, reason, rows) with rows = [(mask, rhs)] affine equations on the
    n * C(n, 2) bits of beta (bit index u((i, j), t) = pair_index * n + t).
    V1 = fixed space of A (the last f coordinates), W1 = fixed space of B (last fB).
      (N1) p odd: V1 is a union of kernel lines, so f is even;
      (N2) p odd, ord_p(2) >= 3: f = fB;
      (N3) p odd: beta(V1, V1) lies in W1 and beta(e, V1) has dimension f - 2 for
           e in V1, so f - 2 <= fB;
      (SB1') basis of V1 in which <e_{n-f+2i}, e_{n-f+2i+1}> are kernel lines;
             for f = 4 the kernel lines in V1 form a regular spread of PG(3, 2),
             normalised to the F_4-spread of omega: e_{n-4} -> e_{n-3} -> e_{n-4}+e_{n-3},
             e_{n-2} -> e_{n-1} -> e_{n-2}+e_{n-1};
      (SB2') beta(e_{n-f}, e_{n-f+j}) = e_{n-fB+j-2} (j = 2..f-1) via GL(W1) in C(B);
             if B = I: beta(e_{n-f}, e_j) = unit vectors for all j outside the line.
      (N5)   f >= 2: y -> beta(e_{n-f}, y) embeds the non-fixed A-module V/V1 into the
             non-fixed B-module (it is injective and intertwines A and B), so the block
             types of A (label la) form a sub-multiset of those of B (label lb);
      (SB2'') and, via the centraliser of B on its non-fixed part, WLOG this embedding
             is the standard one: the first vector of the j-th A-block of type t goes to
             the first vector of the j-th B-block of type t.
    """
    fA2, fB2 = fixed_dim(A, n), fixed_dim(B, n)
    # (N7) any p: for b != 0 annihilating W1, V1 is totally isotropic for B_b (N3), and a
    #      form of rank n-2 has no totally isotropic subspace of dimension > (n+2)/2.
    if fB2 < n and 2 * fA2 > n + 2:
        return False, f"N7: dim V1 = {fA2} > (n+2)/2 while B != I", []
    if p == 2:
        # (N4') B = I: Im(A - I) lies in K_a for every a in V1 \ 0; two different kernel
        #       lines meet in 0, so V1 lies in one kernel line: dim V1 <= 2 < n/2.
        if fB2 == n:
            return False, "N4': B = I with A an involution", []
        # (N10) k = rank(A - I), kB = rank(B - I) (numbers of Jordan blocks of size 2);
        #       U_A = Im(A - I), U_B = Im(B - I). For x in V1 and y = (A - I)y':
        #       beta(x, y) = (B - I) beta(x, y'), so beta(x, U_A) lies in U_B and
        #       k - dim(K_x & U_A) <= kB.
        #       2k < n: x in V1 \ U_A has K_x A-invariant, hence inside V1 (Ay = x + y would
        #       put x in U_A), meeting U_A in dim <= 1: k - 1 <= kB.
        #       2k = n (V1 = U_A): dim(K_x & V1) <= 2 gives kB >= k - 2, with equality only
        #       if every K_x (x in V1) lies in V1, i.e. V1 is a union of kernel lines.
        #       For n = 8 this is kB >= 2; for n = 6 V1 (dim 3) is not a union of lines
        #       (3 does not divide 7), so kB >= k - 1 = 2 as well.
        k, kB = n - fA2, n - fB2
        if 2 * k < n and k - 1 > kB:
            return False, f"N10: involution with k = {k} blocks needs kB >= k - 1 (kB = {kB})", []
        if 2 * k == n and kB < 2:
            return False, f"N10: k = n/2 needs kB >= 2 (kB = {kB})", []
        return True, "", []
    pairs = list(combinations(range(n), 2))
    pidx = {pq: t for t, pq in enumerate(pairs)}

    def u(i, j, t):
        return pidx[(min(i, j), max(i, j))] * n + t

    f, fB = trailing_fixed(A, n), trailing_fixed(B, n)
    if f % 2:
        return False, f"N1: fixed space of A has odd dimension {f}", []
    if ord2(p) >= 3 and f != fB:
        return False, f"N2: dim fixed spaces {f} != {fB}", []
    if f >= 2 and f - 2 > fB:
        return False, f"N3: f - 2 = {f - 2} > fB = {fB}", []
    # (N8) f = 4: V1 is the union of 5 kernel lines (a regular spread of PG(3, 2)), so
    #      beta(x, omega x) = 0 on V1 for its F_4-structure and every component restricted
    #      to V1 has rank 0 or 4; the counting identity for beta restricted to V1 gives
    #      dim span beta(V1, V1) = 2. A B^T-fixed b vanishes on beta(V1, V_nonfixed) and,
    #      if b is orthogonal to beta(V1, V1), on V1 x V1: then R_b contains V1, rank <= n-4.
    #      Such b exist unless fB <= 2; with N3, fB = 2.
    if f == 4 and fB != 2:
        return False, f"N8: f = 4 forces fB = 2 (fB = {fB})", []
    # (N9) p = 3, B fixed-point-free (fB = 0): for x in the omega-part V_w of A,
    #      beta(x, Ax) lies in W1 = 0, so every B_b restricted to V_w is Tr C_b with C_b an
    #      alternating F_4-form on V_w (m_A = dim_{F_4} V_w). If 2 m_A > (n+2)/2, V_w cannot
    #      be totally isotropic for a form of rank n-2, so b -> C_b is injective and
    #      n <= 2 * binom(m_A, 2).
    mA = (n - f) // 2
    if p == 3 and fB == 0 and 4 * mA > n + 2 and n > mA * (mA - 1):
        return False, f"N9: p = 3, B fixed-point-free, m_A = {mA}", []
    # (N6) p = 3, A fixed-point-free (an F_4-structure omega), 4 | n, B != I: for x in V,
    #      beta(x, Ax) = beta(Ax, A^2 x) / B = ... lies in W1, so every b != 0 annihilating
    #      W1 has B_b(x, omega x) = 0 for all x, hence (proof of the Desarguesian theorem)
    #      dim R_b = 0 mod 4, never 2.
    if p == 3 and f == 0 and fB < n and n % 4 == 0:
        return False, "N6: p = 3, A fixed-point-free, B != I, 4 | n", []
    # (SB3') p = 3, f = 0, B = I (n = 8): A-invariant lines are F_4-points of omega = A; the
    #      85 kernel lines are permuted by A in orbits of size 1 or 3 and 85 = 1 mod 3, so
    #      some F_4-point is a kernel line; GL_{F_4}(V) (= C(A)) is transitive on F_4-points:
    #      WLOG <e0, e1> = F_4 e0 is a kernel line; then SB2 through GL(W) = C(I).
    if p == 3 and f == 0 and fB == n and n == 8:
        rows0 = []
        for t in range(n):
            rows0.append((1 << (pidx[(0, 1)] * n + t), 0))
        for j in range(2, n):
            for t in range(n):
                rows0.append((1 << (pidx[(0, j)] * n + t), 1 if t == j - 2 else 0))
        return True, "", rows0
    rows = []

    def vec_eq(terms, target=0):
        """sum of beta(e_i, e_j) over terms (list of (i, j)) equals target vector."""
        for t in range(n):
            mask = 0
            for (i, j) in terms:
                if i != j:
                    mask ^= 1 << u(i, j, t)
            rows.append((mask, target >> t & 1))

    o = n - f
    for i in range(f // 2):
        vec_eq([(o + 2 * i, o + 2 * i + 1)])
    if f == 4:
        def om(k):
            return [k + 1] if (k - o) % 2 == 0 else [k - 1, k]
        for i in range(o, n):
            vec_eq([(i, j) for j in om(i)])
            for k in range(i + 1, n):
                vec_eq([(i, j) for j in om(k)] + [(k, j) for j in om(i)])
    if f >= 2 and la is not None and lb is not None:
        tA = [int(x[1:]) for x in la.split("+")]
        tB = [] if lb == "I" else [int(x[1:]) for x in lb.split("+")]
        for t in set(tA):
            if tA.count(t) > tB.count(t):
                return False, f"N5: A-module {la} does not embed in B-module {lb}", []
        d = ord2(p)
        seen = {}
        for k, t in enumerate(tA):
            j = seen.get(t, 0)
            seen[t] = j + 1
            kb = [i for i, x in enumerate(tB) if x == t][j]
            vec_eq([(o, k * d)], 1 << (kb * d))
    if f >= 2:
        if fB == n:        # B = I: full SB2 through GL(V)
            js = [j for j in range(n) if j not in (o, o + 1)]
            for idx, j in enumerate(js):
                vec_eq([(o, j)], 1 << idx)
        else:
            for j in range(2, f):
                vec_eq([(o, o + j)], 1 << (n - fB + j - 2))
    return True, "", rows


def affine_space(A, B, n, extra=()):
    """Solutions of the equivariance equations plus extra affine rows:
    returns (offset, basis, pairs, u) or None if inconsistent."""
    basis0, pairs, u = equivariant_basis(A, B, n)
    N = len(pairs) * n
    # equivariance rows: recompute as (mask, 0)
    rows = []
    for (i, j) in pairs:
        for t in range(n):
            row = 0
            for (k, l) in pairs:
                coef = ((A[i] >> k & 1) & (A[j] >> l & 1)) ^ ((A[i] >> l & 1) & (A[j] >> k & 1))
                if coef:
                    row ^= 1 << u((k, l), t)
            for s in range(n):
                if B[s] >> t & 1:
                    row ^= 1 << u((i, j), s)
            if row:
                rows.append((row, 0))
    rows += list(extra)
    piv = {}
    for r, rhs in rows:
        for col, (pr, prhs) in piv.items():
            if r >> col & 1:
                r ^= pr
                rhs ^= prhs
        if r:
            col = r.bit_length() - 1
            for c2 in list(piv):
                if piv[c2][0] >> col & 1:
                    piv[c2] = (piv[c2][0] ^ r, piv[c2][1] ^ rhs)
            piv[col] = (r, rhs)
        elif rhs:
            return None
    offset = 0
    for col, (pr, prhs) in piv.items():
        if prhs:
            offset |= 1 << col
    free = [col for col in range(N) if col not in piv]
    basis = []
    for fcol in free:
        v = 1 << fcol
        for col, (pr, prhs) in piv.items():
            if pr >> fcol & 1:
                v |= 1 << col
        basis.append(v)
    return offset, basis, pairs, u


def orbit_reps(M, n, transpose=False):
    """One representative per orbit of <M> (or of <M^T>) on nonzero vectors."""
    def act(v):
        if not transpose:
            return apply(M, v)
        return sum((bin(M[s] & v).count("1") & 1) << s for s in range(n))
    seen, reps = set(), []
    for v in range(1, 1 << n):
        if v in seen:
            continue
        reps.append(v)
        w = v
        while w not in seen:
            seen.add(w)
            w = act(w)
    return reps


def build_sym(n, A, B, perm=True, reduce=False, nf_rows=None):
    """reduce: impose the rank constraints only for one b per <B^T>-orbit and, if
    perm is False, the witness constraints only for one a per <A>-orbit; the others
    follow from beta(Ax, Ay) = B beta(x, y).
    nf_rows: extra affine rows from normal_form() (None: plain equivariant space)."""
    enc = Enc()
    if nf_rows is None:
        basis, pairs, u = equivariant_basis(A, B, n)
        offset = 0
    else:
        sol = affine_space(A, B, n, nf_rows)
        if sol is None:
            raise ValueError("normal-form equations inconsistent")
        offset, basis, pairs, u = sol
    tv = [enc.var(("t", m)) for m in range(len(basis))]
    c = {}
    for pq in pairs:
        c[pq] = [enc.xor([tv[m] for m, v in enumerate(basis) if v >> u(pq, t) & 1]
                         + ([T] if offset >> u(pq, t) & 1 else []))
                 for t in range(n)]
    ell = [[enc.var(("l", i, k)) for k in range(n)] for i in range(n)]
    if reduce:
        b_list = orbit_reps(B, n, transpose=True)
        a_list = None if perm else orbit_reps(A, n)
        add_core(enc, n, c, ell, perm, a_list=a_list, b_list=b_list)
    else:
        add_core(enc, n, c, ell, perm)
    return enc, c, ell, len(basis)


def main():
    cmd, n = sys.argv[1], int(sys.argv[2])
    cs = cases(n)
    if cmd == "list":
        for idx, (p, la, A, lb, B) in enumerate(cs):
            dim = len(equivariant_basis(A, B, n)[0])
            print(f"{idx:3d}  p={p:<4d} A={la:<14s} B={lb:<14s} free beta bits={dim}")
        return
    idx = int(sys.argv[3])
    p, la, A, lb, B = cs[idx]
    perm = "--noperm" not in sys.argv
    t0 = time.time()
    enc, c, ell, dim = build_sym(n, A, B, perm)
    head = (f"n={n} case={idx} p={p} A={la} B={lb} perm={perm} free_beta_bits={dim} "
            f"vars={enc.pool.top} clauses={len(enc.clauses)}")
    if cmd == "cnf":
        with open(sys.argv[4], "w") as f:
            f.write(f"c {head}\n")
            f.write(f"p cnf {enc.pool.top} {len(enc.clauses)}\n")
            for cl in enc.clauses:
                f.write(" ".join(map(str, cl)) + " 0\n")
        print(head, flush=True)
        return
    from pysat.solvers import Solver
    print(head, f"(build {time.time() - t0:.1f}s)", flush=True)
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        t1 = time.time()
        res = s.solve()
        dt = time.time() - t1
        if res:
            C, L = extract(s.get_model(), c, ell, n)
            table = lut(C, L, n)
            ok, du = check(table, n)
            print(f"SAT in {dt:.1f}s: permutation={ok} DU={du}")
            print("LUT:", table)
        else:
            print(f"UNSAT in {dt:.1f}s")


if __name__ == "__main__":
    main()
