/*
 * Size estimates for the enumeration of constant-rank-(n-2) subspaces of Alt(n, F_2).
 * For random constant-rank W_d (built greedily), count exactly the cosets x + W_d all of
 * whose elements have rank n - 2 ("good cosets"), and compare with the independence
 * heuristic 2^(NP-d) p^(2^d), p = fraction of forms of rank n - 2.
 * Usage: crsub_sample n maxd trials [seed] [mincount]  (count only for d >= mincount)
 */
#include "crsub_common.h"
#include <math.h>

int main(int argc, char **argv) {
    init_pairs(atoi(argv[1]));
    int maxd = atoi(argv[2]), trials = atoi(argv[3]);
    int mincount = argc > 5 ? atoi(argv[5]) : 1;
    if (argc > 4) rng_s ^= strtoull(argv[4], 0, 10) * 0x9E3779B97F4A7C15ull;
    long long hist[9];
    build_good(hist);
    printf("n=%d NP=%d ranks:", n, NP);
    for (int r = 0; r <= n; r += 2) printf(" r%d=%lld", r, hist[r]);
    double p = (double)hist[n - 2] / (double)((uint64_t)1 << NP);
    printf("  p=%.5f\n", p);
    static uint32_t el[1 << 8];
    for (int t = 0; t < trials; t++) {
        Sub S; sub_init(&S);
        printf("trial %d:", t);
        for (int d = 1; d <= maxd; d++) {
            /* add a random good coset */
            sub_elements(&S, el);
            int cnt = 1 << S.d, tries = 0;
            uint32_t x;
            do { x = (uint32_t)rnd() & FULL; tries++; } while ((sub_reduce(&S, x) == 0 || !coset_good(el, cnt, x)) && tries < 200000000);
            if (tries >= 200000000) { printf(" [no ext at d=%d]", d); break; }
            sub_add(&S, x);
            sub_elements(&S, el);
            cnt = 1 << S.d;
            if (d < mincount) { printf(" d=%d", d); continue; }
            uint64_t total = (uint64_t)1 << (NP - S.d), goodc = 0;
            for (uint64_t i = 1; i < total; i++) {
                uint32_t rep = sub_unindex(&S, (uint32_t)i);
                if (coset_good(el, cnt, rep)) goodc++;
            }
            double heur = (double)total * pow(p, cnt);
            printf(" d=%d good=%llu (heur %.0f, ratio %.2f)", d, (unsigned long long)goodc, heur, goodc / heur);
            fflush(stdout);
        }
        printf("\n");
    }
    return 0;
}
