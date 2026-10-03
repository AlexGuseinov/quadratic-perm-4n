/*
 * All 6-dimensional subspaces W of Alt(6, F_2) containing J whose nonzero elements all have
 * rank 4 -- WITHOUT any condition on the numbers dim P_a(W) (P_a(W) = elements of W with a in
 * the radical). Same depth-first search as crsub_brute2.c (each W is reached exactly once,
 * through its greedy basis modulo <J>), but a coset x + W_d is accepted as soon as all its
 * elements have rank n - 2; there is no pruning by dim P_a and no full-point condition.
 * For each W found, prints its reduced echelon basis and max_a dim P_a(W); the summary line
 * gives the number of W and the number of W with dim P_a(W) = 2 for all a.
 * Purpose: to decide whether condition (ii) of the paper follows from constant rank for n = 6.
 * Usage: crsub_brute3 n > list.txt        (n = 6; about half an hour)
 */
#include "crsub_common.h"

static uint32_t J, *C1; static size_t nC1; static int Jpiv;
static long long found = 0, nodes[10];

static uint8_t row_of(uint32_t x, int a) {
    uint8_t row = 0;
    for (int p = 0; p < NP; p++) if (x >> p & 1) {
        int i = PI_[p], j = PJ_[p];
        if (a >> i & 1) row ^= 1 << j;
        if (a >> j & 1) row ^= 1 << i;
    }
    return row;
}

/* dim P_a(W) for all a, from all elements (no incremental tricks) */
static void dimp(const uint32_t *el, int cnt, int *dp) {
    for (int a = 1; a < (1 << n); a++) {
        int z = 0;
        for (int t = 1; t < cnt; t++) if (!row_of(el[t], a)) z++;
        dp[a] = z == 0 ? 0 : z == 1 ? 1 : z == 3 ? 2 : z == 7 ? 3 : 4;
    }
}

static void print_canon(const Sub *W) {
    uint32_t b[32]; int d = W->d;
    memcpy(b, W->b, d * 4);
    for (int i = 0; i < d; i++) for (int j = i + 1; j < d; j++) if (b[j] > b[i]) { uint32_t t = b[i]; b[i] = b[j]; b[j] = t; }
    for (int i = 0; i < d; i++) printf("%u%c", b[i], i + 1 < d ? ' ' : '\n');
}

static long long found2 = 0, maxhist[8];

static void dfs(Sub *W, size_t from) {
    nodes[W->d]++;
    uint32_t el[256]; sub_elements(W, el);
    int cnt = 1 << W->d;
    if (W->d == n) {
        int dp[256], mx = 0, all2 = 1;
        dimp(el, cnt, dp);
        for (int a = 1; a < (1 << n); a++) { if (dp[a] > mx) mx = dp[a]; if (dp[a] != 2) all2 = 0; }
        found++; found2 += all2; maxhist[mx]++;
        printf("%d ", mx); print_canon(W); return;
    }
    for (size_t t = from; t < nC1; t++) {
        uint32_t x = C1[t];
        if (!sub_reduce(W, x)) continue;
        /* greedy basis: x must be the smallest representative (mod J) of its coset x + W */
        int ok = 1;
        for (int s = 1; s < cnt && ok; s++) {
            uint32_t y = x ^ el[s];
            if ((y >> Jpiv) & 1) y ^= J;
            if (y < x) ok = 0;
        }
        if (!ok) continue;
        for (int s = 0; s < cnt && ok; s++) if (!good(x ^ el[s])) ok = 0;
        if (!ok) continue;
        Sub W2 = *W; sub_add(&W2, x);
        dfs(&W2, t + 1);
    }
}

int main(int argc, char **argv) {
    init_pairs(atoi(argv[1]));
    long long hist[9]; build_good(hist);
    J = 0;
    for (int k = 0; k < (n - 2) / 2; k++) J |= 1u << PAIR[2 * k][2 * k + 1];
    Sub W; sub_init(&W); sub_add(&W, J);
    Jpiv = W.piv[0];
    C1 = malloc(((size_t)1 << NP) * 4);
    for (uint32_t i = 1; i < (1u << (NP - 1)); i++) {
        uint32_t rep = sub_unindex(&W, i);
        if (good(rep) && good(rep ^ J)) C1[nC1++] = rep;
    }
    for (size_t i = 1; i < nC1; i++) { uint32_t v = C1[i]; size_t j = i; while (j > 0 && C1[j - 1] > v) { C1[j] = C1[j - 1]; j--; } C1[j] = v; }
    setvbuf(stdout, NULL, _IOLBF, 0);
    dfs(&W, 0);
    fprintf(stderr, "constant-rank subspaces containing J: %lld; with dim P_a = 2 for all a: %lld; "
            "by max dim P_a (0..4):", found, found2);
    for (int k = 0; k <= 4; k++) fprintf(stderr, " %lld", maxhist[k]);
    fprintf(stderr, "; nodes by dimension:");
    for (int d = 1; d <= n; d++) fprintf(stderr, " %lld", nodes[d]);
    fprintf(stderr, "\n");
    return 0;
}
