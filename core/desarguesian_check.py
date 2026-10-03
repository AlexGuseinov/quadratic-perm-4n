"""Two checks on sampled quadratic DU-4 permutations of F_2^6 (exact SAT model):

(A) Is the isotropic spread Sigma = {K_a} Desarguesian?  i.e. is there a linear
    omega with omega(x) in K_x \\ {0, x} for all x != 0 (then omega^2+omega+1 = 0
    and Sigma is the set of F_4-points for the F_4-structure given by omega).
    Candidates: omega(e_i) is one of the two other points of K_{e_i}: 2^n choices.

(B) For the Desarguesian ones: Phi(x, y) = beta(omega x, y) + w^2 beta(x, y)
    (values in Y (x) F_4, stored as pairs (u, v) = u + w v) is F_4-bilinear.
"""

import random
from collections import Counter

from pysat.solvers import Solver

from du4perm_sat import build, extract, lut, check


def lin(images, x, n):
    y = 0
    for i in range(n):
        if x >> i & 1:
            y ^= images[i]
    return y


def find_omega(table, n):
    N = 1 << n
    beta = lambda u, v: table[u ^ v] ^ table[u] ^ table[v] ^ table[0]
    K = {a: [v for v in range(N) if beta(a, v) == 0] for a in range(1, N)}
    choices = [[v for v in K[1 << i] if v not in (0, 1 << i)] for i in range(n)]
    found = []
    for m in range(1 << n):
        imgs = [choices[i][(m >> i) & 1] for i in range(n)]
        if all(lin(imgs, x, n) in K[x] and lin(imgs, x, n) not in (0, x)
               for x in range(1, N)):
            found.append(imgs)
    return found, beta


def f4_mul(w, pair):
    """multiply u + w v (w^2 = w + 1) by w: w(u + w v) = w u + (w+1) v = v + w (u + v)"""
    u, v = pair
    return (v, u ^ v)


def check_phi(table, n, imgs, beta):
    N = 1 << n
    om = lambda x: lin(imgs, x, n)

    def phi(x, y):  # beta(om x, y) + w^2 beta(x, y); w^2 = 1 + w
        b = beta(x, y)
        return (beta(om(x), y) ^ b, b)

    for x in range(N):
        for y in range(N):
            if phi(om(x), y) != f4_mul(None, phi(x, y)):
                return False
            if phi(x, om(y)) != f4_mul(None, phi(x, y)):
                return False
    return True


def main():
    n = 6
    enc, c, ell = build(n, sb=True, sb3=True)
    free = [v for p, vs in c.items() for v in vs if isinstance(v, int) and not isinstance(v, bool)]
    stats = Counter()
    with Solver(name="cadical195", bootstrap_with=enc.clauses) as s:
        for t in range(40):
            random.seed(7000 + t)
            assum = [v if random.random() < 0.5 else -v for v in random.sample(free, 8)]
            if not s.solve(assumptions=assum) and not s.solve():
                break
            m = s.get_model()
            C, L = extract(m, c, ell, n)
            tab = lut(C, L, n)
            perm, du = check(tab, n)
            assert perm and du == 4
            oms, beta = find_omega(tab, n)
            desarg = bool(oms)
            phi_ok = desarg and all(check_phi(tab, n, o, beta) for o in oms)
            stats[(desarg, phi_ok, len(oms))] += 1
            pos = set(x for x in m if x > 0)
            s.add_clause([-v if v in pos else v for v in free])
    print("(Desarguesian spread, Phi F4-bilinear for every omega, #omega): count")
    for k, v in sorted(stats.items()):
        print(" ", k, v)


if __name__ == "__main__":
    main()
