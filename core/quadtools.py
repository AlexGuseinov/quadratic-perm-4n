"""Basic tools for quadratic vectorial Boolean functions F: F_2^n -> F_2^n.

Vectors are Python ints (bit i = coordinate i). A function is a lookup table
(list of length 2^n). Everything here is brute force and meant for n <= 8.
"""

from itertools import combinations


def popcount_parity(x: int) -> int:
    return bin(x).count("1") & 1


# ---------------------------------------------------------------- GF(2^n)

def gf_mul(a: int, b: int, n: int, poly: int) -> int:
    """Multiply in GF(2^n) = F_2[X]/(poly); poly includes the X^n term."""
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> n:
            a ^= poly
    return r


def gf_pow(a: int, e: int, n: int, poly: int) -> int:
    r = 1
    while e:
        if e & 1:
            r = gf_mul(r, a, n, poly)
        a = gf_mul(a, a, n, poly)
        e >>= 1
    return r


def power_map(e: int, n: int, poly: int) -> list[int]:
    return [gf_pow(x, e, n, poly) for x in range(1 << n)]


# ---------------------------------------------------------------- linear algebra over F_2

def rank_rows(rows: list[int]) -> int:
    """Rank of a matrix given as a list of row bitmasks."""
    rows = [r for r in rows if r]
    rank = 0
    while rows:
        pivot = rows.pop()
        if not pivot:
            continue
        rank += 1
        low = pivot & -pivot
        rows = [r ^ pivot if r & low else r for r in rows]
        rows = [r for r in rows if r]
    return rank


# ---------------------------------------------------------------- quadratic structure

def is_permutation(F: list[int]) -> bool:
    return len(set(F)) == len(F)


def algebraic_degree_le2(F: list[int], n: int) -> bool:
    """True iff every second-order derivative is constant (degree <= 2)."""
    N = 1 << n
    for a in range(1, N):
        for b in range(a + 1, N):
            c = F[a ^ b] ^ F[a] ^ F[b] ^ F[0]
            for x in range(N):
                if F[x ^ a ^ b] ^ F[x ^ a] ^ F[x ^ b] ^ F[x] != c:
                    return False
    return True


def polar_rows(F: list[int], n: int, b: int) -> list[int]:
    """Rows of the alternating matrix of the polar form of component b.F."""
    rows = []
    for i in range(n):
        ei = 1 << i
        row = 0
        for j in range(n):
            ej = 1 << j
            if i != j and popcount_parity(b & (F[ei ^ ej] ^ F[ei] ^ F[ej] ^ F[0])):
                row |= ej
        rows.append(row)
    return rows


def radical(rows: list[int], n: int) -> list[int]:
    """All a with B(a, .) = 0, i.e. the radical of the form."""
    out = []
    for a in range(1 << n):
        acc = 0
        for i in range(n):
            if (a >> i) & 1:
                acc ^= rows[i]
        if acc == 0:
            out.append(a)
    return out


def ker_sizes(F: list[int], n: int) -> dict[int, int]:
    """|ker L_a| for a != 0, where L_a(x) = F(x+a)+F(x)+F(a)+F(0)."""
    N = 1 << n
    out = {}
    for a in range(1, N):
        c = F[a] ^ F[0]
        out[a] = sum(1 for x in range(N) if F[x ^ a] ^ F[x] ^ c == 0)
    return out


def differential_uniformity(F: list[int], n: int) -> int:
    N = 1 << n
    best = 0
    for a in range(1, N):
        counts = [0] * N
        for x in range(N):
            counts[F[x ^ a] ^ F[x]] += 1
        best = max(best, max(counts))
    return best


def sign(perm: list[int]) -> int:
    """+1 for even, -1 for odd permutation (via cycle count)."""
    N = len(perm)
    seen = [False] * N
    cycles = 0
    for s in range(N):
        if not seen[s]:
            cycles += 1
            x = s
            while not seen[x]:
                seen[x] = True
                x = perm[x]
    return 1 if (N - cycles) % 2 == 0 else -1


def component_profile(F: list[int], n: int) -> dict:
    """Ranks and radicals of all nonzero components."""
    ranks = {}
    rads = {}
    for b in range(1, 1 << n):
        rows = polar_rows(F, n, b)
        ranks[b] = rank_rows(rows)
        rads[b] = frozenset(radical(rows, n))
    return {"ranks": ranks, "radicals": rads}


# ---------------------------------------------------------------- alternating forms as bitmasks

def alt_index(n: int) -> list[tuple[int, int]]:
    """Coordinate order for Alt(n, F_2): pairs (i, j), i < j."""
    return list(combinations(range(n), 2))


def alt_rows(mask: int, n: int, idx: list[tuple[int, int]]) -> list[int]:
    rows = [0] * n
    for k, (i, j) in enumerate(idx):
        if (mask >> k) & 1:
            rows[i] |= 1 << j
            rows[j] |= 1 << i
    return rows


def alt_rank(mask: int, n: int, idx: list[tuple[int, int]]) -> int:
    return rank_rows(alt_rows(mask, n, idx))


def count_alt_by_rank_formula(n: int, r2: int, q: int = 2) -> int:
    """Number of n x n alternating matrices over F_q of rank r2 (= 2r)."""
    r = r2 // 2
    num = 1
    for i in range(2 * r):
        num *= q ** (n - i) - 1
    den = 1
    for i in range(r):
        den *= q ** (2 * i + 2) - 1
    val = num
    for i in range(r):
        val *= q ** (2 * i)
    assert val % den == 0
    return val // den
