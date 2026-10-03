"""Sanity checks for the Track D framework (quadratic DU-4 permutations).

1. Count of alternating matrices by rank: formula vs brute force (n = 4, 6),
   then the n = 8 numbers used in the heuristic.
2. x^5 on GF(64): the known n = 6 quadratic DU-4 permutation.
   Checks: permutation, degree 2, DU, component ranks, |ker L_a|,
   radicals forming a line spread, sign.
3. n = 4: maximal dimension of a subspace of Alt(4, F_2) whose nonzero
   elements all have rank 2 (the reason no quadratic DU-4 permutation
   exists at n = 4).
"""

from collections import Counter
from quadtools import (
    alt_index, alt_rank, count_alt_by_rank_formula, power_map, is_permutation,
    algebraic_degree_le2, differential_uniformity, component_profile,
    ker_sizes, sign,
)


def section(title):
    print("\n" + "=" * 70 + "\n" + title + "\n" + "=" * 70)


# ---------------------------------------------------------------- 1
section("1. Alternating matrices over F_2 by rank")
for n in (4, 6):
    idx = alt_index(n)
    brute = Counter(alt_rank(m, n, idx) for m in range(1 << len(idx)))
    formula = {r: count_alt_by_rank_formula(n, r) if r else 1
               for r in range(0, n + 1, 2)}
    print(f"n={n}: brute={dict(sorted(brute.items()))}")
    print(f"      formula={formula}  match={dict(brute) == formula}")
n = 8
f8 = {r: count_alt_by_rank_formula(n, r) if r else 1 for r in range(0, 9, 2)}
total = sum(f8.values())
print(f"n=8: formula={f8}  total={total} (2^28={1 << 28})")
p6 = f8[6] / (1 << 28)
print(f"      fraction of rank-6 forms = {p6:.4f}")

# ---------------------------------------------------------------- 2
section("2. x^5 on GF(2^6), modulus X^6+X+1")
n, poly = 6, 0b1000011
F = power_map(5, n, poly)
print("permutation:", is_permutation(F))
print("degree <= 2:", algebraic_degree_le2(F, n))
print("DU:", differential_uniformity(F, n))
prof = component_profile(F, n)
print("component ranks:", dict(Counter(prof["ranks"].values())))
ks = ker_sizes(F, n)
print("|ker L_a| distribution:", dict(Counter(ks.values())))
rads = prof["radicals"]
distinct = set(rads.values())
print("distinct radicals:", len(distinct),
      " sizes:", dict(Counter(len(r) for r in distinct)))
# spread check: pairwise intersections trivial and union covers F_2^n
pts = Counter(p for r in distinct for p in r if p)
disjoint = all(len(r1 & r2) == 1 for r1 in distinct for r2 in distinct if r1 != r2)
print("radicals pairwise meet only in 0:", disjoint,
      " every nonzero point covered once:", set(pts.values()) == {1} and len(pts) == 63)
# components sharing a radical: do they form a 2-dim space minus 0?
by_rad = {}
for b, r in rads.items():
    by_rad.setdefault(r, []).append(b)
ok = all(len(bs) == 3 and bs[0] ^ bs[1] == bs[2] for bs in
         (sorted(v) for v in by_rad.values()))
print("each radical shared by exactly 3 components {b1, b2, b1+b2}:", ok)
print("sign:", "even" if sign(F) == 1 else "odd")

# ---------------------------------------------------------------- 3
section("3. n=4: max dim of constant-rank-2 subspaces of Alt(4,F_2)")
n = 4
idx = alt_index(n)
rank2 = {m for m in range(1, 1 << len(idx)) if alt_rank(m, n, idx) == 2}
print("rank-2 forms:", len(rank2))

best = 0


def extend(span, elems):
    """span: set of the subspace's elements; try to add a new generator."""
    global best
    dim = len(span).bit_length() - 1
    best = max(best, dim)
    for g in sorted(rank2):
        if g in span or g <= elems:
            continue
        new = {g ^ w for w in span}
        if new <= rank2:
            extend(span | new, g)


extend({0}, 0)
print("maximal dimension:", best, "(needed for a DU-4 permutation: 4)")
