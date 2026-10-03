/*
 * Common helpers for constant-rank subspaces of Alt(n, F_2) (n <= 8).
 *
 * A form is a mask of NP = n(n-1)/2 bits; bit p(i,j) (i < j, pairs in lexicographic order,
 * the same order as itertools.combinations(range(n), 2)) is B(e_i, e_j).
 * GOOD is a bitset over all 2^NP forms: bit set iff rank B = n - 2.
 * Group elements g in GL(n, 2) are stored as n row bytes (row i = linear functional
 * x -> (g x)_i).  The action is B^g(x, y) = B(g x, g y), a right action: B^{gh} = (B^g)^h,
 * where gh is the matrix product.
 */
#ifndef CRSUB_COMMON_H
#define CRSUB_COMMON_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>

static int n, NP, PAIR[8][8], PI_[28], PJ_[28];
static uint64_t *GOOD;          /* bitset, 2^NP bits */
static uint32_t FULL;           /* (1 << NP) - 1 */

static inline int good(uint32_t x) { return (GOOD[x >> 6] >> (x & 63)) & 1; }

static int rank_mask(uint32_t m) {
    uint8_t r[9] = {0};
    for (int p = 0; p < NP; p++)
        if (m >> p & 1) { r[PI_[p]] |= 1 << PJ_[p]; r[PJ_[p]] |= 1 << PI_[p]; }
    int rank = 0;
    for (int col = 0; col < n; col++) {
        int piv = -1;
        for (int i = rank; i < n; i++) if (r[i] >> col & 1) { piv = i; break; }
        if (piv < 0) continue;
        uint8_t t = r[piv]; r[piv] = r[rank]; r[rank] = t;
        for (int i = 0; i < n; i++) if (i != rank && (r[i] >> col & 1)) r[i] ^= r[rank];
        rank++;
    }
    return rank;
}

static void init_pairs(int nn) {
    n = nn; NP = n * (n - 1) / 2; FULL = (NP == 32) ? 0xffffffffu : ((1u << NP) - 1);
    int p = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) { PAIR[i][j] = PAIR[j][i] = p; PI_[p] = i; PJ_[p] = j; p++; }
}

/* rank distribution in hist[0..n]; builds GOOD */
static void build_good(long long hist[9]) {
    size_t words = ((size_t)1 << NP) / 64 + 1;
    GOOD = calloc(words, 8);
    memset(hist, 0, 9 * sizeof(long long));
    for (uint64_t m = 0; m <= FULL; m++) {
        int r = rank_mask((uint32_t)m);
        hist[r]++;
        if (r == n - 2) GOOD[m >> 6] |= 1ull << (m & 63);
    }
}

/* GOOD from a cache file if present (else build and save); hist only filled when built */
static void build_good_cached(const char *path, long long hist[9]) {
    size_t words = ((size_t)1 << NP) / 64 + 1;
    FILE *f = path ? fopen(path, "rb") : NULL;
    if (f) {
        GOOD = malloc(words * 8);
        size_t got = fread(GOOD, 8, words, f);
        fclose(f);
        if (got == words) { memset(hist, 0, 9 * sizeof(long long)); return; }
        free(GOOD);
    }
    build_good(hist);
    if (path && (f = fopen(path, "wb"))) { fwrite(GOOD, 8, words, f); fclose(f); }
}

/* ---------- group elements ---------- */
typedef struct { uint8_t r[8]; } Mat;

static Mat mat_id(void) { Mat g; memset(&g, 0, sizeof g); for (int i = 0; i < n; i++) g.r[i] = 1 << i; return g; }
static Mat mat_mul(Mat a, Mat b) {           /* (ab)x = a(b x): row i of ab = sum_k a_ik row_k(b) */
    Mat c; memset(&c, 0, sizeof c);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++) if (a.r[i] >> k & 1) c.r[i] ^= b.r[k];
    return c;
}
static int mat_inv(Mat a, Mat *out) {
    Mat b = mat_id();
    for (int col = 0; col < n; col++) {
        int piv = -1;
        for (int i = col; i < n; i++) if (a.r[i] >> col & 1) { piv = i; break; }
        if (piv < 0) return 0;
        uint8_t t = a.r[piv]; a.r[piv] = a.r[col]; a.r[col] = t;
        t = b.r[piv]; b.r[piv] = b.r[col]; b.r[col] = t;
        for (int i = 0; i < n; i++) if (i != col && (a.r[i] >> col & 1)) { a.r[i] ^= a.r[col]; b.r[i] ^= b.r[col]; }
    }
    *out = b; return 1;
}
static int mat_eq(Mat a, Mat b) { return memcmp(a.r, b.r, 8) == 0; }
static int mat_is_id(Mat a) { Mat e = mat_id(); return mat_eq(a, e); }
static uint64_t mat_key(Mat a) { uint64_t k; memcpy(&k, a.r, 8); return k; }

static uint32_t wedge(uint8_t u, uint8_t v) {
    uint32_t m = 0;
    for (int p = 0; p < NP; p++) {
        int k = PI_[p], l = PJ_[p];
        if (((u >> k & 1) & (v >> l & 1)) ^ ((u >> l & 1) & (v >> k & 1))) m |= 1u << p;
    }
    return m;
}

/* action tables: 4 chunks of 7 bits */
typedef struct { uint32_t t[4][128]; } Act;
static void act_build(Mat g, Act *A) {
    uint32_t img[28];
    for (int p = 0; p < NP; p++) img[p] = wedge(g.r[PI_[p]], g.r[PJ_[p]]);
    for (int c = 0; c < 4; c++) {
        A->t[c][0] = 0;
        for (int v = 1; v < 128; v++) {
            int low = __builtin_ctz(v), p = 7 * c + low;
            A->t[c][v] = A->t[c][v & (v - 1)] ^ (p < NP ? img[p] : 0);
        }
    }
}
static inline uint32_t act(const Act *A, uint32_t x) {
    return A->t[0][x & 127] ^ A->t[1][(x >> 7) & 127] ^ A->t[2][(x >> 14) & 127] ^ A->t[3][(x >> 21) & 127];
}
static uint32_t act_slow(Mat g, uint32_t x) {
    uint32_t y = 0;
    for (int p = 0; p < NP; p++) if (x >> p & 1) y ^= wedge(g.r[PI_[p]], g.r[PJ_[p]]);
    return y;
}

/* ---------- subspaces of Alt: fully reduced echelon basis (pivot = highest bit) ---------- */
typedef struct { int d; uint32_t b[32]; int piv[32]; uint32_t pivmask; } Sub;
static void sub_init(Sub *S) { S->d = 0; S->pivmask = 0; }
static uint32_t sub_reduce(const Sub *S, uint32_t x) {
    for (int k = 0; k < S->d; k++) if (x >> S->piv[k] & 1) x ^= S->b[k];
    return x;
}
static int sub_add(Sub *S, uint32_t x) {    /* returns 0 if x in S */
    x = sub_reduce(S, x);
    if (!x) return 0;
    int p = 31 - __builtin_clz(x);
    for (int k = 0; k < S->d; k++) if (S->b[k] >> p & 1) S->b[k] ^= x;
    S->b[S->d] = x; S->piv[S->d] = p; S->d++; S->pivmask |= 1u << p;
    return 1;
}
static inline uint32_t sub_index(const Sub *S, uint32_t canon) { return _pext_u32(canon, FULL & ~S->pivmask); }
static inline uint32_t sub_unindex(const Sub *S, uint32_t idx) { return _pdep_u32(idx, FULL & ~S->pivmask); }
/* all 2^d elements of S */
static void sub_elements(const Sub *S, uint32_t *out) {
    out[0] = 0;
    for (int k = 0; k < S->d; k++)
        for (int t = 0; t < (1 << k); t++) out[(1 << k) + t] = out[t] ^ S->b[k];
}
/* coset x + S entirely rank n-2 ? */
static int coset_good(const uint32_t *el, int cnt, uint32_t x) {
    for (int t = 0; t < cnt; t++) if (!good(x ^ el[t])) return 0;
    return 1;
}

/* xorshift RNG */
static uint64_t rng_s = 88172645463325252ull;
static inline uint64_t rnd(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }

#endif
