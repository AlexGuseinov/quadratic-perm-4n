"""MacWilliams test for constant-rank subspaces of Alt(n, F_q).

Setting: Alt(n, F_q) with the inner product <A, B> = sum_{i<j} a_ij b_ij and the
additive character (-1)^<A,B> (q = 2).  Classes of the alternating-forms
association scheme: rank 2i, i = 0..m, m = floor(n/2).

For an additive subgroup W with rank distribution (A_0, ..., A_m), the dual
W^perp has rank distribution
      A'_k = (1/|W|) * sum_i A_i * Q_k(i),
      Q_k(i) = sum_{B : rank B = 2k} (-1)^<A,B>   (any A of rank 2i).
A'_k must be non-negative integers.  A constant-rank-r subspace of dimension d
has A_0 = 1, A_{r/2} = q^d - 1, others 0.

Q_k(i) is computed two ways:
  * brute force (n <= 6, and n = 8 through the C helper eigen_alt8),
  * the Delsarte–Goethals generalized Krawtchouk formula.
"""

import sys
from fractions import Fraction
from itertools import combinations


def qbinom(a, b, q):
    if b < 0 or b > a:
        return 0
    num = den = 1
    for i in range(b):
        num *= q ** (a - i) - 1
        den *= q ** (i + 1) - 1
    return num // den


def dg_eigen(n, q=2):
    """Delsarte–Goethals eigenvalues Q_k(i) of the alternating-forms scheme."""
    m = n // 2
    b = q * q
    c = q ** (n * (n - 1) // (2 * m)) if m else 1
    Q = [[0] * (m + 1) for _ in range(m + 1)]
    for k in range(m + 1):
        for i in range(m + 1):
            s = 0
            for j in range(k + 1):
                s += ((-1) ** (k - j) * b ** ((k - j) * (k - j - 1) // 2)
                      * qbinom(m - j, m - k, b) * qbinom(m - i, j, b) * c ** j)
            Q[k][i] = s
    return Q


def rank_rows(rows):
    rows = [x for x in rows if x]
    rk = 0
    while rows:
        p = rows.pop()
        rk += 1
        low = p & -p
        rows = [x ^ p if x & low else x for x in rows]
        rows = [x for x in rows if x]
    return rk


def brute_eigen(n):
    pairs = list(combinations(range(n), 2))
    P = len(pairs)

    def rows_of(mask):
        rows = [0] * n
        for k, (i, j) in enumerate(pairs):
            if mask >> k & 1:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
        return rows

    ranks = [rank_rows(rows_of(mk)) // 2 for mk in range(1 << P)]
    m = n // 2
    reps = {}
    for mk in range(1 << P):
        reps.setdefault(ranks[mk], mk)
    Q = [[0] * (m + 1) for _ in range(m + 1)]
    for i, A in reps.items():
        for B in range(1 << P):
            Q[ranks[B]][i] += -1 if bin(A & B).count("1") & 1 else 1
    return Q


def dual_distribution(n, r, d, Q, q=2):
    m = n // 2
    A = [0] * (m + 1)
    A[0] = 1
    A[r // 2] = q ** d - 1
    return [Fraction(sum(A[i] * Q[k][i] for i in range(m + 1)), q ** d)
            for k in range(m + 1)]


def verdict(dist):
    ok = all(x.denominator == 1 and x >= 0 for x in dist)
    return "consistent" if ok else "IMPOSSIBLE"


def main():
    print("Check: Delsarte–Goethals formula vs brute force")
    for n in (3, 4, 5, 6):
        bf, dg = brute_eigen(n), dg_eigen(n)
        print(f"  n={n}: match={bf == dg}")
        if bf != dg:
            print("   brute:", bf, "\n   formula:", dg)
            sys.exit(1)

    print("\nConstant rank n-2, dimension d:  dual rank distribution A'_0..A'_m")
    for n in range(4, 21, 2):
        Q = dg_eigen(n)
        for d in range(n - 2, n + 3):
            dist = dual_distribution(n, n - 2, d, Q)
            shown = [str(x) for x in dist]
            print(f"  n={n:2d} (n mod 4 = {n % 4}) d={d:2d}: {verdict(dist):10s} {shown}")


if __name__ == "__main__":
    main()
