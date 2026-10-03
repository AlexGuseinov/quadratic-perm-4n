/*
 * Second check of crsub_enum at n = 6 that does not use the radical-closure rule.
 * Enumerates all 6-dim W containing J with (i) all nonzero elements of rank 4 and
 * (ii) dim P_a(W) = 2 for all a, by a plain depth-first search over admissible cosets:
 * W_{d+1} = W_d + <x> with x taken in increasing order among the representatives modulo <J>
 * (the greedy minimal basis of W / <J> is increasing, so every W is reached), with
 * x also the smallest representative modulo <J> of its coset x + W_d, so that every W is reached
 * exactly once (through its greedy basis), where
 * admissible = all elements of rank n-2 and no element has a full point of W_d in its
 * radical (a point a is full if dim P_a(W_d) = 2; such a point cannot lie in the radical of an
 * element of W outside W_d, by (ii)). Also prunes if dim P_a(W_d) > 2 for some a.
 * Prints the reduced echelon basis of each W found (pipe through sort -u) and counts.
 * Usage: crsub_brute2 n > list.txt
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

static void dfs(Sub *W, size_t from) {
    nodes[W->d]++;
    uint32_t el[256]; sub_elements(W, el);
    int cnt = 1 << W->d, dp[256];
    dimp(el, cnt, dp);
    for (int a = 1; a < (1 << n); a++) if (dp[a] > 2) return;
    if (W->d == n) {
        for (int a = 1; a < (1 << n); a++) if (dp[a] != 2) return;
        found++; print_canon(W); return;
    }
    int full[256], nf = 0;
    for (int a = 1; a < (1 << n); a++) if (dp[a] == 2) full[nf++] = a;
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
        for (int s = 0; s < cnt && ok; s++) {
            uint32_t y = x ^ el[s];
            if (!good(y)) ok = 0;
            for (int k = 0; k < nf && ok; k++) if (!row_of(y, full[k])) ok = 0;
        }
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
    fprintf(stderr, "leaves passing (i),(ii): %lld; nodes by dimension:", found);
    for (int d = 1; d <= n; d++) fprintf(stderr, " %lld", nodes[d]);
    fprintf(stderr, "\n");
    return 0;
}
