# Working lemmas (Mövzu 2, Track D) — status: proved here, not yet checked by Ali

**Attribution (added 2026-09-28).** L1 is Charpin–Peng, FFA 56 (2019) 188–208, Theorem 5
(written there as sum_lambda 2^{l(lambda)} = sum_a 2^{d(a)}). L2 is Charpin–Peng Theorem 1 +
Lemma 4 + Theorem 7 (quadratic DU-4 permutations are two-valued {0,4}, no bent components,
NL = 2^{n-1} - 2^{n/2}, single amplitude). Charpin–Peng Theorem 8: for DO polynomials with all
exponents 2^i + 2^j, t | i, j, t | n, and no bent component, two-valued {0, 2^t} iff every
radical has dimension t (they note radicals are F_{2^t}-subspaces) — no condition on n/t.
Our new part: L3/L4 (kernel spread), L6 parity (m odd), L7, L8.

Setting: F: V = F_2^n -> F_2^n quadratic, F(0) = 0, bilinear (polar) map
beta(u, v) = F(u+v) + F(u) + F(v). Component b: B_b(u, v) = <b, beta(u, v)>,
an alternating form with radical R_b. For a != 0: ker_a = {v : beta(a, v) = 0},
s_a = dim ker_a (>= 1 since a in ker_a).

## L1 (count identity)
sum_{a != 0} (2^{s_a} - 1) = sum_{b != 0} (2^{n - rank B_b} - 1).
Both sides count pairs (a, b), a != 0, b != 0, with B_b(a, .) = 0.

## L2 (DU-4 permutation => constant rank n-2 and s_a = 2)
If F is a permutation, every component is balanced, hence not bent, so
rank B_b <= n-2 and the right side of L1 is >= 3(2^n - 1).
DU(F) = max_a 2^{s_a}; DU = 4 gives the left side <= 3(2^n - 1).
So equality: rank B_b = n-2 for all b and s_a = 2 for all a.
(For non-permutations this fails: x^5 on GF(2^8) has DU 4 with bent and
rank-4 components.)

## L3 (isotropic spread)
If s_a = 2 for all a != 0, the sets ker_a are 2-dim subspaces forming a
line spread Sigma of V, and beta vanishes on L x L for every L in Sigma.
Proof: u in ker_a \ {0, a} => a in ker_u (beta symmetric), so
ker_u contains <a, u> = ker_a, equal by dimension.
Checked on 30 random n = 6 SAT solutions (logs/structure_n6.log).
Note: the radicals R_b need NOT be the spread lines (true for x^5 on
GF(2^6), false for 29 of 30 random n = 6 solutions).

## L4 (spread to spread)
Under L3, F is linear on every L in Sigma (beta = 0 there, F(0) = 0).
If F is also a permutation, F(L) is a 2-dim subspace and {F(L)} is a line
spread Sigma' of the codomain.

## L5 (n = 4, MacWilliams)
In Alt(4, F_2) no 4-dim subspace has all nonzero elements of rank 2: the
MacWilliams transform in the alternating-forms scheme gives a negative
entry (-2) (logs/macwilliams_alt.log). For n = 8 and d = 8 the transform is
non-negative and integral, so MacWilliams alone does not settle n = 8.

## L6 (F_4-quadratic family: no DU-4 permutation when 4 | n) — proof by hand
Let n = 2m and F: F_{2^n} -> F_{2^n} be quadratic with F_4-bilinear beta
(exactly the maps sum c_ij x^{4^i + 4^j} + l(x), l any F_2-linear map). Every F_2-component is Tr_{F_4/F_2}(g_b) with
g_b = Tr_{F_{2^n}/F_4}(b F), an F_4-quadratic form whose polar form B_g is
F_4-bilinear and alternating. The F_2-radical of Tr(g_b) equals the
F_4-radical of B_g (Tr_{F_4/F_2} is onto and B_g(x, .) is F_4-linear). An
alternating form over F_4 has even rank, so the F_4-radical has dimension
m - (even), and the F_2-radical dimension is 2(m - even). For m even this
is a multiple of 4, never 2. So no component has rank n - 2, and by L2 no
member of this family is a DU-4 permutation when n = 0 mod 4.
For m odd, radical dimension 2 is possible (Gold x^5 on GF(2^6) is an example).
Check: 40 random members at n = 8 give only ranks {0, 4, 8}
(logs/f4_homogeneous_check.log).
Status: short and elementary; likely known in the "APN over F_q" literature —
must be checked before claiming novelty. Same argument over F_{2^t}: in the
F_{2^t}-quadratic class, two-valued permutations with delta = 2^t need n/t odd.

## L7 (Desarguesian isotropic spread is impossible when 4 | n) — proof by hand
Suppose there is a linear omega on V with omega^2 + omega + 1 = 0 and
beta(x, omega x) = 0 for all x (equivalently: the isotropic spread of L3 is
the Desarguesian spread of the F_4-structure given by omega).
(i) Polarising, beta(x, omega y) = beta(omega x, y): omega is self-adjoint.
(ii) Put Phi(x, y) = beta(omega x, y) + w^2 beta(x, y) in Y (x) F_4 (w in F_4,
     w^2 = w + 1). Using (i) and omega^2 = omega + 1 one checks
     Phi(omega x, y) = w Phi(x, y) = Phi(x, omega y); Phi is symmetric and
     Phi(x, x) = 0. So Phi is F_4-bilinear alternating on V (F_4-dim m = n/2).
(iii) beta is an F_2-linear image of Phi (its w-coordinate), so every component
     B_b = Tr_{F_4/F_2} o C with C an F_4-bilinear alternating form, and as in L6
     rad B_b = rad C has F_2-dimension 2(m - rank_{F_4} C), a multiple of 4 when m
     is even. So no component has rank n - 2, contradicting L2.
Consequence: for n = 0 mod 4, the isotropic spread of a quadratic DU-4
permutation is non-Desarguesian. L6 is the special case beta F_4-bilinear.
n = 4: every line spread of PG(3,2) is Desarguesian, which re-proves D(4).
Checks: x^5 on GF(2^6) has two such omega and Phi is F_4-bilinear
(desarguesian_check.py). SAT for n = 8 with the Desarguesian constraint:
UNSAT in 1 s, DRAT proof VERIFIED by drat-trim (logs/drat_d8_desarg.log).
Observation: 40 of 40 sampled DU-4 permutations at n = 6 have
NON-Desarguesian isotropic spreads (logs/desarguesian_check_n6.log). This is
n = 2 mod 4 data; it does not tell us whether non-Desarguesian kernel spreads
can occur when 4 | n (that is exactly the open case).

## L8 (general t — unifies L7 with the classical APN fact) — proof by hand
Let F be a quadratic permutation of F_2^n with s_a = t for all a != 0 (i.e. every
derivative kernel K_a has dimension t; DU = 2^t, two-valued). Suppose there is an
F_{2^t}-structure on V (V = F_{2^t}^m, n = tm) with K_x = F_{2^t} x for all x != 0.
[CORRECTED after independent review: for t >= 3 the K_a need NOT form a spread
 from s_a = t alone (the L3 argument uses t = 2); without the permutation
 hypothesis this is false (reviewer's examples n=6, t=3 and t=4). So the spread
 structure is now part of the hypothesis.]
Claim: m = n/t is odd.
Proof: beta(x, lambda x) = 0 for all lambda in F_{2^t}; polarising, every lambda
acts self-adjointly: beta(lambda x, y) = beta(x, lambda y). Descent: with
{theta_j} a basis of F_{2^t} over F_2 and {theta_j*} its trace-dual basis,
Phi(x, y) = sum_j beta(theta_j x, y) (x) theta_j*  (values in Y (x) F_{2^t})
is F_{2^t}-bilinear, symmetric, alternating, and beta = (Tr (x) id) o Phi.
[For t = 2 this is the explicit Phi of L7.] Hence each component
B_b = Tr o C_b with C_b F_{2^t}-bilinear alternating, and rad B_b = rad C_b has
F_2-dimension t(m - 2 r_b) (C_b has even F_{2^t}-rank 2 r_b).
If m is even, every radical dimension is a multiple of 2t (0 = bent is excluded
because F is a permutation, so components are balanced), so
sum_b (2^{dim rad} - 1) >= (2^n - 1)(2^{2t} - 1) > (2^n - 1)(2^t - 1) =
sum_a (2^{s_a} - 1), contradicting L1. So m is odd.
* t = 1: the spread is the set of points, always "Desarguesian" over F_2: this is
  the classical fact that quadratic APN permutations need n odd.
* t = 2: L7 (and L6 as a special case).
* Sharp: Gold x^{2^k+1} with gcd(k, n) = t and n/t odd.
Descent step written in full in paper/main.tex (Theorem general-t); checked
numerically by the reviewer for (n,t) = (6,2), (8,2), (6,3), (9,3), (8,4).

## Observations (n = 6 samples; not proofs)
* The Pluecker vectors u^v of the isotropic spread lines span a space of
  dimension C(n,2) - n = 9 in all 30 samples, i.e. Sigma^perp = W exactly:
  the spread determines the component space (logs/spread_span.log).
  (MacWilliams forces W^perp to contain exactly (2^n-1)/3 rank-2 forms, i.e.
  the common totally isotropic lines of W are exactly the spread lines.)
  Desarguesian spread: span 9 (n = 6), 16 (n = 8, Sigma^perp of dim 12).
* Self-adjoint algebra {phi : phi^T M = M phi for all M in W}: dimension 1
  (scalars only) in all 30 samples, dimension 2 (= F_4) for x^5
  (logs/selfadjoint_n6.log). So the non-Desarguesian examples carry no hidden
  global F_4-structure.

* Sub-spreads: number of pairs of spread lines whose 4-dim span contains 5
  spread lines — 210 of 210 for x^5 (Desarguesian), exactly 50 in all 30
  samples (logs/spread_closed_pairs_n6.log). Suggests one non-Desarguesian
  spread type at n = 6 (with 5 sub-spreads of PG(3,2)). To be identified.
* Local version of L7: on a closed 4-space U (5 spread lines, necessarily a
  Desarguesian spread of U) the Phi-argument applies to beta|_{U x U}, so every
  component restricted to U has rank 0 or 4, never 2. (Proof sketch; check.)

* Single-form test (negative): a line spread totally isotropic for ONE
  alternating form of rank 2r exists for (n, 2r) = (4,4), (6,4), (6,6), (8,6);
  not for (4,2) (logs/ti_spread_small.log, logs/ti_spread_n8_r6.log). So the
  obstruction, if any, needs the whole 8-dim W, not a single component.

## SB1 (symmetry breaking on the input side) — WLOG for DU-4 permutations
The isotropic spread Sigma of L3 contains n/2 lines L_1, ..., L_{n/2} with
V = L_1 + ... + L_{n/2} (direct sum). For n = 8: two distinct spread lines
are skew, U = L_1 + L_2 has dim 4. If no spread line were skew to U, the
240 points outside U would all lie on lines meeting U in one point, i.e.
on at most 15 lines carrying 2 outside points each (30 < 240), a
contradiction; so some L_3 is skew to U. The same count with
U_3 = U + L_3 (63 points; 192 points outside, at most 63 lines with 2
outside points each = 126 < 192) gives L_4. General n: same argument,
since 2 * (2^k - 1) < 2^n - 2^k for every k <= n - 2.
Choosing the basis e_{2i}, e_{2i+1} of L_{i+1} gives beta(e_{2i}, e_{2i+1}) = 0.

## SB2 (symmetry breaking on the output side) — WLOG
ker beta(e_0, .) = L_1 = <e_0, e_1> (L3), so beta(e_0, e_j), j = 2..n-1, are
linearly independent. An output change of basis (A o F, A in GL(n, 2)),
which preserves DU and bijectivity, maps them to the first n-2 unit vectors.
SB1 and SB2 act on different sides and are compatible.

## Analogy (motivation for the conjecture)
t = 1 (APN): quadratic APN permutations exist iff n is odd.
t = 2 (DU 4, two-valued): Gold x^{2^k+1} with gcd(k, n) = 2 is a permutation
iff n/2 is odd. Conjecture D: quadratic DU-4 permutations exist iff n = 2 mod 4.

## Independent review (2026-09-27, separate agent, re-derived all identities)
Verdicts: L1, L2, L4, SB1/SB2, L6, L7 correct; L3 correct with two added sentences
(partition, beta(L,L)=0 — added to the paper); L8 conclusion and descent correct,
but the t-spread claim for t >= 3 was unproved/false without permutation -> removed
from the statement (now a hypothesis). Novelty: kernel isotropy is in MTX19 Lemma 1
(cited); MTX19 Thm 5 gives delta >= q for F_q-quadratic maps but no parity / m odd
statement (checked). n = 4: every line spread of PG(3,2) is Desarguesian (reviewer
enumerated all 56); needs a citation. Wording fixes applied in main.tex: n = 8 DRAT
certificate is only a consistency check of the encoding (the case is proved by L7);
omega-normalisation stated; n = 6 remark no longer claims anything about 4 | n.

## Symmetric case (2026-09-28) — proofs by hand, used by sym_search.normal_form
Setting: F quadratic, s_a = 2 for all a != 0 and all components of rank n-2 (true for
DU-4 permutations). Suppose beta(Ax, Ay) = B beta(x, y) with A of odd prime order p and
B^p = I. V1 = ker(A - I), f = dim V1; W1 = ker(B - I), fB = dim W1 (= dim ker(B^T - I)).
(S0) K_{Aa} = A K_a;  Im beta(Aa, .) = B Im beta(a, .);  R_{B^T b} = A^{-1} R_b.
     (from beta(Aa, Av) = B beta(a, v) and B_b(Ax, Ay) = B_{B^T b}(x, y).)
(N1) a in V1 \ 0 => K_a subset of V1.  Proof: A K_a = K_a, K_a = {0, a, y, a+y};
     Ay in {y, a+y}; Ay = a+y gives A^2 y = y, and with A^p y = y (p odd) Ay = y,
     contradiction. So V1 is a union of kernel lines and f is even.
(N3) x, y in V1 => beta(x, y) = beta(Ax, Ay) = B beta(x, y), so beta(V1, V1) subset W1;
     for a in V1 \ 0, beta(a, .) restricted to V1 has kernel K_a (dim 2), so its image
     has dim f - 2 <= fB.
(N2) If ord_p(2) >= 3 (p >= 5), then f = fB.  Proof: nontrivial irreducible F_2<A>-
     modules have dimension ord_p(2) >= 3, so every A-invariant (B^T-invariant) 2-dim
     subspace lies in the fixed space. For a in V1, Im beta(a, .)^perp = {b : a in R_b}
     is B^T-invariant of dim 2, so all 3 such b are B^T-fixed; for B^T-fixed b, R_b is
     A-invariant, so R_b subset V1. Count pairs (a, b), a in V1 \ 0, b B^T-fixed, a in R_b:
     3(2^f - 1) = 3(2^fB - 1).
(N4) If B = I then f = 0.  Proof: for a in V1 \ 0, beta(a, (A - I)y) = (B - I) beta(a, y)
     = 0, so Im(A - I) subset K_a subset V1 (N1); but Im(A - I) meets V1 only in 0
     (A semisimple, p odd), so A = I, contradiction.
Normalisations (WLOG, changes of basis in C(A) and C(B)):
(SB1') GL(V1) commutes with A: choose a basis of V1 made of f/2 kernel lines in direct
     sum (the SB1 counting argument inside V1); for f = 4 the 5 kernel lines in V1 form a
     spread of PG(3, 2), which is regular; its stabiliser acts transitively on ordered
     pairs of its lines, so WLOG it is the F_4-spread of omega with e_{n-3} = omega e_{n-4},
     e_{n-1} = omega e_{n-2}, i.e. beta(x, omega x) = 0 on V1.
(SB2') GL(W1) commutes with B: the f - 2 independent vectors beta(e_{n-f}, e_{n-f+j})
     (j = 2..f-1) lie in W1 (N3) and are mapped to the first f - 2 basis vectors of W1.
     If B = I: all six beta(e_{n-f}, e_j), j outside the kernel line, to unit vectors.
Validation: at n = 6 every case gives the same SAT/UNSAT verdict with and without these
conditions (51 cases; script in the session scratchpad, results in logs/sym_n6*.jsonl).
Consequence for n = 8 (final, with N5-N10 below): 79 of the 93 pairs are excluded by the
lemmas; the remaining 14 are settled by computation (logs/sym_summary_n8.txt).
(N5)  f >= 2: for x in V1, y -> beta(x, y) on V_nf = Im(A - I) is injective (kernel K_x
      lies in V1) and intertwines A and B, so V_nf embeds in W_nf = Im(B - I) as a module;
      WLOG (centraliser of B on W_nf) the embedding for x = e_{n-f} is the standard one.
(N6)  p = 3, f = 0, B != I, 4 | n: for x in V, beta(x, Ax) = beta(Ax, A^2 x) (A^2 = A + 1)
      = B beta(x, Ax), so beta(x, Ax) lies in W1; for b != 0 orthogonal to W1,
      B_b(x, omega x) = 0 for all x, and the proof of Theorem desarg gives dim R_b = 0 mod 4.
(N7)  any p, B != I: V1 is totally isotropic for B_b, b != 0 orthogonal to W1 (N3); a
      form of rank n-2 has no totally isotropic subspace of dim > (n+2)/2: f <= (n+2)/2.
(N8)  f = 4 (odd p): V1 is a union of 5 kernel lines (regular spread of PG(3,2)); local
      counting identity for beta on V1 x V1 (all components of rank 0 or 4 there) gives
      dim span beta(V1, V1) = 2; a B^T-fixed b orthogonal to it has R_b containing V1:
      so fB <= 2, and fB = 2 by N3.
(N9)  p = 3, fB = 0: beta(x, Ax) = 0 on the omega-part V_w of A (dim_{F4} = m_A), so
      B_b = Tr C_b on V_w with C_b in Alt(m_A, F_4); if 4 m_A > n + 2 then C_b != 0 for
      b != 0 (N7-type bound) and n <= 2 binom(m_A, 2).
(N4') p = 2, B = I: Im(A - I) lies in every K_x, x in V1 \ 0: impossible for n >= 6.
(N10) p = 2, k = rank(A - I), kB = rank(B - I): for x in V1, beta(x, U_A) lies in U_B
      (U = Im(M - I)); x in V1 \ U_A has K_x inside V1 meeting U_A in dim <= 1, so
      k - 1 <= kB when 2k < n; when 2k = n (V1 = U_A): kB >= k - 2, with equality only if
      V1 is a union of kernel lines (n = 8: kB >= 2; kB = 2 forces the union).
(N8 detail) the counting identity is used for the quadratic map V1 -> V (codomain V):
      15 (2^{8-2} - 1) = 15 #{b != 0 : B_b vanishes on V1 x V1}, so beta(V1, V1) spans 2 dims.
Split cases (scripts sym_case44.py, sym_case32.py, sym_p2.py; justifications in their
docstrings): p = 5 fixed-point-free (8 alternatives for K_{e0}), p = 3 case 32 (42
conjugacy classes of M = phi_{e6}^{-1} phi_{e7} in GL(3, 4)), p = 2 cases (2 alternatives
for the kernel line through a chosen point; F_4-normalised V1 when kB = 2).

## General n = 8 case: local structure of the kernel spread (2026-09-28)
Setting: beta with s_a = 2 for all a != 0 and all components of rank n - 2 (every a lies in
exactly 3 radicals R_b, namely b in Im beta(a, .)^perp). For kernel lines L, M let
U = L + M (dim 4) and k(L, M) = number of kernel lines inside U.
(K4) k != 4 (any n). [Simpler proof, used in the paper: sum_{x in U} q(x) = 0 for
     q(x) = beta(x, omega x), since each component is a quadratic form on F_2^4 and takes
     the value 1 an even number of times; so q_L + q_M + q_N + q_P + q_Q = 0.] The 5 lines of the regular spread of PG(U) through L, M (and a third
     kernel line N if k >= 3) are the F_4-points of some omega on U. For each b,
     Q_b(x) = B_b(x, omega x) is a quadratic form on U, constant on each F_4-point
     (Q_b(omega x) = B_b(omega x, x + omega x) = Q_b(x)). It vanishes on the kernel lines in
     U. If L, M, N are kernel lines, the zeros are 10 + 3[q_P = 0] + 3[q_Q = 0] for the other
     two F_4-points P, Q; a quadratic form on F_2^4 has 0, 4, 6, 8, 10, 12 or 16 zeros, so
     q_P = q_Q for every b: P is a kernel line iff Q is (a line is a kernel line iff
     beta(x, omega x) = 0 on it). Hence k in {2, 3, 5}; for k = 3, beta(p, omega p) = v is
     the same nonzero vector for all 6 points p of P and Q.
(K5) No closed pair (k = 5) when n = 8. If k = 5, omega is self-adjoint for beta on U, so
     every B_b restricted to U has rank 0 or 4. The counting identity on U
     (Remark rem:count-general) gives 15 * 63 = 15 * #{b != 0 : B_b|U = 0}, so
     S = span beta(U, U) has dimension 2. For each of the 63 values b != 0 orthogonal to S,
     U is totally isotropic for B_b; a form of rank 6 on F_2^8 has maximal totally
     isotropic subspaces of dimension 5 = 2 + 3, so U meets R_b nontrivially. Counting
     pairs (a, b) with a in U \ 0, b in S^perp \ 0, a in R_b: at least 63, at most
     3 * 15 = 45. Contradiction. (For n = 6 the same count gives equality: R_b lies in U
     for the 15 such b; closed pairs do occur at n = 6. Checked on 450 closed pairs of
     sampled n = 6 solutions: check_closed_lemma.py.)
     In particular this re-proves Theorem desarg for n = 8 (there every pair is closed).
(K6) General bound from the same argument: for every subspace W of dim w, each b != 0 with
     B_b|W = 0 has dim(W & R_b) >= w - 3 (n = 8); so no 6-dim W is totally isotropic for
     any component, i.e. span beta(W, W) = F_2^8 for every 6-dim W.
Remaining cases (n8_cases.py): B1 (some pair with k = 3) and B2 (all pairs k = 2).
At n = 6: B1 is UNSAT (also without the permutation condition), B2 is SAT.
(K7) Pair types (n = 8): every pair with k = 2 has d = dim span beta(L, M) = 4, every pair with
     k = 3 has d = 3. Proof: per b, counting pairs (x, y) with B_b(x, y) = 0 gives
     3 n0(b) + n1(b) = 5292 (n_r = ordered pairs of kernel lines on which B_b has rank r);
     summing over b: sum over unordered pairs of (3 N0(U) + N2(U)) = 189 * 3570. The local
     identity gives 3N0 + N2 = 157 for (k,d) = (2,3), 189 for (2,4) and (3,3), and (3,4) is
     impossible (N2 <= 128 when k = 3). So every term equals 189. Checked at n = 6 (the
     analogous identity gives d = 4 for all k = 2 pairs): 1760 pairs (k=2, d=4), 550 (k=5, d=2).
(K7 general) The pair-type identity holds for every even n >= 6 (reviewer's remark, now in
     the paper as Prop. prop:local (b)): 3 n0(b) + n1(b) = (2^n - 4)(2^{n-1} - 2)/6 per b,
     average of 3N0 + N2 over pairs = 3 * 2^{n-2} - 3; k = 2 => d = 4, k = 3 => d = 3,
     k = 5 => d = 2.
Status of the general n = 8 case (2026-09-28, midday):
  * closed pairs impossible (proved); k in {2, 3}; d = 4 for k = 2, d = 3 for k = 3.
  * B1 (some k = 3 pair): explicit normal form of beta on U = L + M (N = graph of id):
    beta(e0,e2) = e'0, beta(e0,e3) = e'1, beta(e1,e2) = e'1, beta(e1,e3) = e'0 + e'1 + v,
    q = v on P u Q; restrictions of components to U = symmetric 2x2 matrices in these bases
    (classes: 31 rank 0 [b _|_ S], 96 F4-type rank 4, 96 rank 2 with radical one of the three
    transversals, 32 rank 4). Counting consequences found so far are all consistent
    (e.g. number of points of L u M u N with v in Im_a equals 2 + #(alpha-radicals inside U)).
    At n = 6 B1 is excluded by the radical count (R_b inside U for the 7 b _|_ S: 21 < 27).
  * B2 (all pairs k = 2): no contradiction found yet.
  * SAT: n8_cases.py (B1, B2, --local, --d4); runs in logs/kissat_d8_B*.log.

## n = 8 by exhaustive search over component spaces (2026-09-28, afternoon)
Reduction. For a quadratic DU-4 permutation F of F_2^n, W = {B_b : b in F_2^n} is an n-dim
subspace of Alt(n, F_2) with
  (i)  every nonzero element of rank n - 2 (L2), and
  (ii) dim P_a = 2 for every a != 0, where P_a = {w in W : a in rad w} (P_a = Im beta(a,.)^perp,
       dim = s_a = 2).
More generally, (i)+(ii) hold for every quadratic F whose derivative kernels all have dim 2 and
which has no bent component (by L1 the ranks are then all n - 2). crsub_enum.c searches all
such W for n = 8 up to GL(8, 2).
Radical-closure lemma (completeness). Let W be a solution containing a subspace W_d and let r be
a point with dim(P_r cap W_d) <= 1. Then P_r (dim 2 in W) is not inside W_d, so W contains an
element c outside W_d with r in rad c. The coset c + W_d lies in W, so all its elements have rank
n - 2, and no element of it has in its radical a point a that is "full" for W_d
(dim P_a cap W_d = 2), since otherwise dim P_a >= 3. So c + W_d is a candidate in
Cand(W_d, r) := {cosets x + W_d : all elements rank n-2, no full point of W_d in any radical,
some element has r in its radical}, and W contains W_d + <c>.
Search. W_1 = <J> (J = e0^e1 + e2^e3 + e4^e5; all rank-6 forms are equivalent). At W_d choose r
(d <= 2: first point with dim P_r cap W_d = 1, else first with 0; d >= 3: the point minimising
4 * |Cand| + [dim = 0]); children W_d + c for c in Cand(W_d, r), one per orbit of a group that
stabilises W_d and fixes r. Pruning (all valid by the lemma and counting inside W / W_d):
  * dim(P_a cap W_d) <= 2 for all a;
  * |C(W_d)| >= 2^(n-d) - 1 (C = good cosets without full points in radicals);
  * every a with dim(P_a cap W_d) = 1 needs >= 1 candidate coset, with dim 0 needs >= 3.
Symmetry: G_1 = Stab(J) from explicit generators (checked); stabilisers of points and cosets
from random products and random Schreier generators (a subgroup is enough for correctness);
every group element used is checked at run time to stabilise W_d and fix r.
Validation.
  * n = 6: the search finds exactly 9 distinct W (27 leaves), of two GL-types: the Desarguesian
    one (x^5; 210 closed pairs) and the non-Desarguesian one with 50 closed pairs (the type of
    all 30 sampled SAT solutions); both carry DU-4 permutations (stage 2). Without symmetry
    (--nosym): 12288 distinct W containing J, all of these two types (invariants checked on a
    sample of 1500).
  * n = 8, independent implementation (Python + numpy, written by a separate agent from the
    mathematical specification only, without seeing the C code): rank-6 table recomputed by
    Pfaffians, identical to the C table (149,920,960 forms); node counts of the level-3
    subtrees 0, 5, 17 (without symmetry) identical to the C program at every level, including
    the pruning statistics; the certificate for the symmetry reductions at levels 1-2
    (cert8.txt, 1,468,672 lines: for every candidate coset an explicit matrix mapping an orbit
    representative onto it, stabilising W_d and fixing r) checked completely; the 486 roots
    W_3 reproduced as a set of subspaces.
First full run (2026-09-28, 12:58-15:05, 2 processes, logs/crsub_n8_part{0,1}.log): all 486
level-3 nodes processed (hashes agree with the separate listing), no BUG line, empty solution
files. Nodes: dimension 3: 486, dimension 4: 1,466,318, dimension 5: 244,686,904, dimension 6: 0.
Every node of dimension 5 fails (P2) (244,668,896) or (P3) (18,008). About 4.1 CPU-hours.
=> NO SOLUTION: no quadratic DU-4 permutation of F_2^8.
Independent referee (separate agent, 2026-09-28, /tmp review_n8.md -> summarised here): no error in
the mathematics or the code; asked for (B1) a rerun from archived source with provenance
(run_n8.sh, logs/crsub_n8_final_*; sources copied to logs/crsub_n8_final_src/), (B2) a log line
for every level-3 index (added: "pruned" lines), (B3) the independent Python program on all 486
roots (it searches without symmetry below level 3, so the level-3 orbit reduction of the C
program needs no certificate), (B4) cross-checks under all three level-2 nodes (roots 243 and 480
added: C --nosym-from 3 and Python agree exactly), text fixes (B5), and the rank-table hash (B6,
in the provenance file).
n = 6 consistency (referee): the stabilisers of the two classes have orders 362,880 (= |GammaL(3,4)|,
x^5 space) and 5,760; |GL(6,2)| * 63 / (|Stab| * 18,228) gives 192 and 12,096 subspaces
containing J, total 12,288 = the no-symmetry count. crsub_brute2.c (rule-free DFS over admissible
cosets, each subspace exactly once through its greedy basis mod J; 21 min): exactly the same
12,288 subspaces (logs/crsub_brute2_n6.log, logs/crsub_brute2_n6_list.txt).
Final runs (2026-09-28 evening).
  * Rerun from the archived sources (run_n8.sh; provenance in logs/crsub_n8_final_provenance.txt):
    identical numbers at every dimension; check_n8_logs.py: PASS (logs/crsub_n8_final_check.txt).
  * Independent Python implementation on all 486 roots, no symmetry below dimension 3
    (second_impl/batch_results_*.txt, batch_summary.txt): 8,604,140 nodes of dim 4,
    1,381,518,812 of dim 5 (1,381,416,881 fail P2, 101,931 fail P3), none of dim 6, 0 solutions,
    13.4 CPU-hours. Together with the certificate of levels 1-2 this is a second complete proof.
    (8,604,140 = the number of candidates of the 486 level-3 nodes in the C program, 8.604e6.)
Paper, Appendix B (app:search, added 2026-09-28 late evening): the search as a procedure
Search(U, H) (steps 1-5), Lemma B-action (group action preserves nodes, admissible cosets and
candidates), Lemma B-step (incremental formulas used by the programs: dim P_a(U') = dim P_a(U) +
eps_a; admissible cosets of U' from those of U; candidate test; existence of a point with
dim <= 1), Theorem B-complete (completeness, by induction on d), implementation notes with the
correspondence to crsub_enum.c functions, and the table of node counts of both programs.
Paper revision 2026-09-29: Section 6 shortened to Proposition prop:local (with proof) and a
summary of the earlier approaches (structured classes, SAT, symmetry theorem); the removed text
and the old Appendix A are in paper/removed_2026-09-29_section6_appendixA.tex. The search
appendix is now Appendix A (label app:search). Added citation De Meyer-Bilgin, ToSC 2019(2)
169-192 (classification of quadratic permutations for n <= 6; eight classes with delta = 4 for
n = 6 -- number taken from their Table 15 as reported by a web fetch; to be checked by the author).
Removed the uncited reference C17. AI-use statement drafted (author to adapt).
Paper revision 2026-10-02 (after three referee-style reports and one more internal review):
  * Section 4 rebuilt around Theorem general-t with Lemma trace-rad (rad(Tr o C) = rad C);
    Desarguesian obstruction and the F_4 statement are corollaries; Corollary cor:general-t
    reduced to its two statements, commentary moved to a remark.
  * Remark after Theorem constrank: the claimed equivalence "single amplitude <=> two-valued" for
    quadratic functions without bent components was wrong as stated (the direction single
    amplitude => two-valued is the open question of Section 7); now: two-valued => single
    amplitude, and the converse if delta <= 4, both by the counting identity.
  * Section 5: "point" defined (nonzero vector), letters for forms and points fixed, (P1) no
    longer listed as a pruning test (it holds for every node), CPU times stated as measured
    (two processes on a two-core VM), second implementation described as separately written
    (not "independent"); directory renamed independent/ -> second_impl/.
  * Appendix A rewritten as "The computation in mathematical terms": objects (F_2^28 model,
    coordinate order, rank table by two methods, coset representatives, phi_a), the two lemmas,
    procedure and completeness theorem, how each step is computed, symmetry and certificate
    (full checklist), what the proof depends on (table), n = 6 worked example, why the search
    is small (heuristic vs measured counts, logs/crsub_sample_n8.log), counts.
  * New computation, n = 6: crsub_brute3.c enumerates ALL 6-dim constant-rank-4 subspaces
    containing J (no condition on the P_a): 12,288 subspaces, identical to the list of
    crsub_brute2, all with dim P_a = 2 (logs/crsub_brute3_n6.log; nodes by dimension
    1, 5086, 1344240, 14232960, 3251712, 12288). So for n = 6 condition (ii) follows from (i),
    and Alt(6, F_2) has exactly two GL(6,2)-classes of such subspaces. A second, different
    check (reduction to a 3-dim subspace with a common radical point, then exhaustive
    extension) gave the same answer: 0 of 15,600 such 3-dim subspaces extend.
    The earlier sentence "does not follow from the constant-rank condition alone" was
    unsupported and is replaced.
  * The certificate of levels 1-2 is now shipped (second_impl/data/cert8.txt.xz, 9.8 MB);
    verify_cert.py reads the compressed file.
  * "Use of AI tools" section removed from the paper (the author states this to the journal).
  * Open for the author before submission: check the citations DMB19 ("eight classes with
    delta = 4" for n = 6), MTX19 Lemma 1 and Theorem 5 against the published versions; make the
    repository public (the paper gives its URL) and archive a release.
  * Possible extension: Conjecture H for n = 8 (no 8-dim constant-rank-6 subspace at all). For
    d <= 6 a point with dim P_r(U) = 0 always exists (3(2^d - 1) < 255), so the radical-point
    prescription works without assuming (ii); only admissibility by full points is lost.
    Rough size estimate 10^9 - 10^10 nodes: tens to hundreds of CPU-hours.
Fourth referee-style report (2026-10-02, afternoon): no mathematical objection; "major revision"
only because the repository URL gives 404. Changes made:
  * Artifact: ARTIFACT.md (trusted base, commands, expected output, determinism),
    expected_counts.json (compared automatically by check_n8_logs.py and
    second_impl/summarize_batch.py), second_impl/run_all.sh (one command for the verification
    that does not use crsub_enum; steps 1-2 and one subtree re-tested from a clean data
    directory: PASS, root 0 with the same counts). task1_pfaffian.py no longer needs good8.bin;
    verify_cert.py and summarize_batch.py return an exit status.
  * Paper: threshold-implementation motivation made precise (d+1 shares, uniformity separate);
    "quadratic" = degree at most 2; Desarguesian t-spread defined in Theorem general-t (t | n);
    "partial normal form"; explicit paragraph after Theorem n8 and in the conclusion that the
    theorem does not prove Conjecture H; isotropic-dimension bound 2 + 3 explained; run times
    given per process (2.2 + 2.1 h; rerun 3.4 + 3.3 h on the shared machine; second program
    6.8 + 6.7 h); determinism of the random choices stated; Reproducibility section restructured.
  * New references: Gupta-Mandal, arXiv:2310.03340 (Galois-theoretic constant-rank subspaces;
    characteristic not 2 throughout, so not applicable) and Qu-Xiong-Li, FFA 24 (2013) 55-65
    (negative answer to Bracken-Tan-Tan's Problem 5.1: differentially 4-uniform but not always
    permutations).
  * Still to do by the author: make the repository public and deposit a release with a DOI
    (e.g. Zenodo); then replace the URL sentence in the Reproducibility section by the DOI.

Folder reorganisation (2026-10-03). Earlier entries of this file name files as they were in the
flat layout. New places: crsub_*.c, crsub_common.h, crsub_invariants.py, run_n8.sh,
check_n8_logs.py -> search_n8/; quadtools.py, du4perm_sat.py, constrank_sat.py,
validate_basics.py, macwilliams_alt.py, check_x5.py, desarguesian_check.py -> core/;
*_n6.py, regulus_switch.py, check_closed_lemma.py -> dimension6/; constrank8.c ->
earlier_approaches/structured/; export_dimacs.py, du4perm_cms.py, n8_cases.py, cube_*.py,
run_d8_half.py, ti_spread*.py, vps/, cnf/, proofs/ -> earlier_approaches/sat/; sym_*.py,
beta_enum.c, validate_sym_n6.py, conjH_sym.py -> earlier_approaches/symmetry/; bin/,
build_mac.sh -> earlier_approaches/; NOTES_lemmas.md, REVIEW_*.md -> notes/.
Checked after the move: sha256 of crsub_enum.c, crsub_common.h and of the rebuilt binary equal
the provenance record; check_all.sh; crsub_enum 6 + crsub_invariants.py (2 classes);
check_n8_logs.py and summarize_batch.py PASS; sym_summary.py reproduces the table of
logs/sym_summary_n8.txt; one short run of each script in sat/, symmetry/ and dimension6/.
