"""Hand-checkable examples for L3 and L7.

1. x^5 on GF(16): kernels K_a = {x : beta(a, x) = 0} equal a*F_4 (L3 example;
   not a permutation, L3 does not need one).
2. x^5 on GF(64): permutation with DU 4; kernels form a Desarguesian spread;
   the map Phi of L7 is F_4-bilinear.
"""

from quadtools import gf_mul, gf_pow, power_map
from desarguesian_check import find_omega, check_phi
from du4perm_sat import check


def example_gf16():
    n, poly = 4, 0b10011  # X^4 + X + 1
    F = power_map(5, n, poly)
    beta = lambda u, v: F[u ^ v] ^ F[u] ^ F[v] ^ F[0]
    F4 = [x for x in range(1 << n) if gf_pow(x, 3, n, poly) == 1] + [0]  # 0,1,w,w^2
    ok = True
    for a in range(1, 1 << n):
        K = sorted(v for v in range(1 << n) if beta(a, v) == 0)
        aF4 = sorted(gf_mul(a, c, n, poly) for c in F4)
        ok &= (K == aF4)
    perm, du = check(F, n)
    print(f"[1] x^5 on GF(16): K_a == a*F_4 for all a: {ok};  permutation={perm}, DU={du}")


def example_gf64():
    n, poly = 6, 0b1000011  # X^6 + X + 1
    F = power_map(5, n, poly)
    perm, du = check(F, n)
    oms, beta = find_omega(F, n)
    phi_ok = all(check_phi(F, n, o, beta) for o in oms)
    print(f"[2] x^5 on GF(64): permutation={perm}, DU={du}; "
          f"Desarguesian omegas found={len(oms)}; Phi F_4-bilinear={phi_ok}")


if __name__ == "__main__":
    example_gf16()
    example_gf64()
