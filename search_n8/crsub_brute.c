/*
 * Independent check of crsub_enum for small n (n = 6): all n-dim subspaces W of Alt(n, F_2)
 * containing J, all nonzero elements of rank n-2, with dim P_a = 2 for every a != 0.
 * Plain DFS over good cosets with increasing representatives (no radical rule, no groups):
 * the greedy minimal basis of W / <J> is increasing, so every W is reached; duplicates are
 * removed by printing the reduced echelon basis (pipe through sort -u).
 * Usage: crsub_brute n > list.txt
 */
#include "crsub_common.h"

static uint32_t J, *C1; static size_t nC1;
static long long found = 0;

static int check_ii(const Sub *W) {
    uint32_t el[256]; sub_elements(W, el);
    int cnt = 1 << W->d;
    for (int a = 1; a < (1 << n); a++) {
        int z = 0;
        for (int t = 1; t < cnt; t++) {
            uint32_t x = el[t]; int row = 0;
            for (int p = 0; p < NP; p++) if (x >> p & 1) {
                int i = PI_[p], j = PJ_[p];
                if (a >> i & 1) row ^= 1 << j;
                if (a >> j & 1) row ^= 1 << i;
            }
            if (!row) z++;
        }
        if (z != 3) return 0;
    }
    return 1;
}

static void print_canon(const Sub *W) {
    uint32_t b[32]; int d = W->d;
    memcpy(b, W->b, d * 4);
    for (int i = 0; i < d; i++) for (int j = i + 1; j < d; j++) if (b[j] > b[i]) { uint32_t t = b[i]; b[i] = b[j]; b[j] = t; }
    for (int i = 0; i < d; i++) printf("%u%c", b[i], i + 1 < d ? ' ' : '\n');
}

/* every point lies in the radicals of at most 3 nonzero elements (dim P_a <= 2) */
static int dimp_ok(const uint32_t *el, int cnt) {
    for (int a = 1; a < (1 << n); a++) {
        int z = 0;
        for (int t = 1; t < cnt; t++) {
            uint32_t x = el[t]; int row = 0;
            for (int p = 0; p < NP; p++) if (x >> p & 1) {
                int i = PI_[p], j = PJ_[p];
                if (a >> i & 1) row ^= 1 << j;
                if (a >> j & 1) row ^= 1 << i;
            }
            if (!row && ++z > 3) return 0;
        }
    }
    return 1;
}

static void dfs(Sub *W, size_t from) {
    if (W->d == n) { if (check_ii(W)) { found++; print_canon(W); } return; }
    uint32_t el[256]; sub_elements(W, el);
    int cnt = 1 << W->d;
    if (W->d >= 3 && !dimp_ok(el, cnt)) return;
    for (size_t t = from; t < nC1; t++) {
        uint32_t x = C1[t];
        if (!sub_reduce(W, x)) continue;
        if (!coset_good(el, cnt, x)) continue;
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
    C1 = malloc(((size_t)1 << NP) * 4);
    for (uint32_t i = 1; i < (1u << (NP - 1)); i++) {
        uint32_t rep = sub_unindex(&W, i);
        if (good(rep) && good(rep ^ J)) C1[nC1++] = rep;
    }
    /* C1 is increasing in the compressed index; sort by value for the greedy-basis argument */
    for (size_t i = 1; i < nC1; i++) { uint32_t v = C1[i]; size_t j = i; while (j > 0 && C1[j - 1] > v) { C1[j] = C1[j - 1]; j--; } C1[j] = v; }
    dfs(&W, 0);
    fprintf(stderr, "leaves passing (i),(ii): %lld\n", found);
    return 0;
}
