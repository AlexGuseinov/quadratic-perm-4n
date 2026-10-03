/*
 * Constant-rank search in Alt(8, F_2) for quadratic DO maps over GF(2^8).
 *
 * Question tested (hypothesis H, n = 8): is there a quadratic F: F_2^8 -> F_2^8
 * all of whose 255 nonzero components have rank exactly 6?  This is necessary
 * for a quadratic permutation with differential uniformity 4.
 *
 * F(x) = sum_{i<j} c_ij x^(2^i + 2^j) over GF(2^8), modulus X^8+X^4+X^3+X^2+1.
 * Component b: f_b(x) = Tr(b F(x)); its polar form is the 8x8 alternating
 * matrix M_b[u][v] = Tr(b * beta(e_u, e_v)), beta the symmetric bilinear map
 * of F.  "bad" = number of b != 0 with rank(M_b) != 6.
 *
 * Modes:
 *   validate                  known power maps
 *   binomial                  x^d1 + c x^d2, all d1<d2, c != 0
 *   trinomial [maxbad]        x^d1 + c2 x^d2 + c3 x^d3 (exhaustive)
 *   f2coef [maxbad]           all 2^28 F_2-coefficient DO polynomials
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 8
#define NPAIR 28
#define POLY 0x11D

static uint8_t *RANK;            /* rank of alternating matrix, 2^28 entries */
static uint8_t MUL[256][256];
static uint8_t TR[256];
static int PI[NPAIR], PJ[NPAIR]; /* pair index -> (i, j), i < j */
static uint8_t T[NPAIR][NPAIR];  /* T[term ij][pair uv] = beta_ij(e_u, e_v) */
static uint32_t MC[NPAIR][256][N]; /* component masks of c*x^(d_t) */

static uint8_t gmul(uint8_t a, uint8_t b) {
    uint16_t r = 0, aa = a;
    while (b) {
        if (b & 1) r ^= aa;
        b >>= 1;
        aa <<= 1;
        if (aa & 0x100) aa ^= POLY;
    }
    return (uint8_t)r;
}

static uint8_t gpow(uint8_t a, int e) {
    uint8_t r = 1;
    while (e) {
        if (e & 1) r = MUL[r][a];
        a = MUL[a][a];
        e >>= 1;
    }
    return r;
}

static int rank_rows(uint8_t rows[N]) {
    uint8_t r[N];
    memcpy(r, rows, N);
    int rank = 0;
    for (int col = 0; col < N; col++) {
        int piv = -1;
        for (int k = rank; k < N; k++)
            if (r[k] >> col & 1) { piv = k; break; }
        if (piv < 0) continue;
        uint8_t t = r[piv]; r[piv] = r[rank]; r[rank] = t;
        for (int k = 0; k < N; k++)
            if (k != rank && (r[k] >> col & 1)) r[k] ^= r[rank];
        rank++;
    }
    return rank;
}

static void mask_to_rows(uint32_t m, uint8_t rows[N]) {
    memset(rows, 0, N);
    for (int k = 0; k < NPAIR; k++)
        if (m >> k & 1) {
            rows[PI[k]] |= 1u << PJ[k];
            rows[PJ[k]] |= 1u << PI[k];
        }
}

static void init(void) {
    int k = 0;
    for (int i = 0; i < N; i++)
        for (int j = i + 1; j < N; j++) { PI[k] = i; PJ[k] = j; k++; }
    for (int a = 0; a < 256; a++)
        for (int b = 0; b < 256; b++) MUL[a][b] = gmul(a, b);
    for (int a = 0; a < 256; a++) {
        uint8_t x = a, t = 0;
        for (int i = 0; i < N; i++) { t ^= x; x = MUL[x][x]; }
        TR[a] = t; /* t is 0 or 1 */
    }
    /* beta_ij(x, y) = x^(2^i) y^(2^j) + x^(2^j) y^(2^i) */
    for (int t = 0; t < NPAIR; t++)
        for (int p = 0; p < NPAIR; p++) {
            uint8_t eu = 1u << PI[p], ev = 1u << PJ[p];
            int i = PI[t], j = PJ[t];
            T[t][p] = MUL[gpow(eu, 1 << i)][gpow(ev, 1 << j)] ^
                      MUL[gpow(eu, 1 << j)][gpow(ev, 1 << i)];
        }
    RANK = malloc((size_t)1 << NPAIR);
    uint8_t rows[N];
    for (uint32_t m = 0; m < (1u << NPAIR); m++) {
        mask_to_rows(m, rows);
        RANK[m] = (uint8_t)rank_rows(rows);
    }
}

/* component masks for output basis b_k = x^k, from beta(e_u,e_v) values */
static void comp_masks(const uint8_t beta[NPAIR], uint32_t masks[N]) {
    for (int k = 0; k < N; k++) {
        uint32_t m = 0;
        uint8_t b = 1u << k;
        for (int p = 0; p < NPAIR; p++)
            if (TR[MUL[b][beta[p]]]) m |= 1u << p;
        masks[k] = m;
    }
}

/* number of nonzero components whose rank != 6; stops early above maxbad */
static int count_bad(const uint32_t masks[N], int maxbad, int hist[N + 1]) {
    uint32_t cur = 0;
    int bad = 0;
    for (uint32_t g = 1; g < 256; g++) {
        int bit = __builtin_ctz(g);
        cur ^= masks[bit];
        int r = RANK[cur];
        if (hist) hist[r]++;
        if (r != 6 && ++bad > maxbad && !hist) return bad;
    }
    return bad;
}

static void build_mc(void) {
    for (int t = 0; t < NPAIR; t++)
        for (int c = 0; c < 256; c++) {
            uint8_t beta[NPAIR];
            for (int p = 0; p < NPAIR; p++) beta[p] = MUL[c][T[t][p]];
            comp_masks(beta, MC[t][c]);
        }
}

static void beta_from_terms(const int *terms, const uint8_t *coef, int nt,
                            uint8_t beta[NPAIR]) {
    memset(beta, 0, NPAIR);
    for (int q = 0; q < nt; q++)
        for (int p = 0; p < NPAIR; p++)
            beta[p] ^= MUL[coef[q]][T[terms[q]][p]];
}

static int exponent(int t) { return (1 << PI[t]) + (1 << PJ[t]); }

static void report_power(int d) {
    /* locate term index for x^d */
    for (int t = 0; t < NPAIR; t++)
        if (exponent(t) == d) {
            uint8_t beta[NPAIR], one = 1;
            uint32_t masks[N];
            int hist[N + 1] = {0};
            beta_from_terms(&t, &one, 1, beta);
            comp_masks(beta, masks);
            count_bad(masks, 1000, hist);
            printf("x^%-3d ranks: r0=%d r2=%d r4=%d r6=%d r8=%d\n", d,
                   hist[0], hist[2], hist[4], hist[6], hist[8]);
        }
}

int main(int argc, char **argv) {
    const char *mode = argc > 1 ? argv[1] : "validate";
    int maxbad = argc > 2 ? atoi(argv[2]) : 3;
    init();
    build_mc();
    long cnt[256] = {0}; /* histogram of exact counts <= maxbad */
    if (!strcmp(mode, "validate")) {
        long hist_all[N + 1] = {0};
        for (uint32_t m = 0; m < (1u << NPAIR); m++) hist_all[RANK[m]]++;
        printf("Alt(8,F2) by rank: r0=%ld r2=%ld r4=%ld r6=%ld r8=%ld\n",
               hist_all[0], hist_all[2], hist_all[4], hist_all[6], hist_all[8]);
        int ds[] = {3, 5, 9, 17, 33, 65, 129, 6, 10, 12};
        for (unsigned q = 0; q < sizeof ds / sizeof *ds; q++) report_power(ds[q]);
        return 0;
    }
    if (!strcmp(mode, "binomial")) {
        int best = 1 << 30;
        long hits = 0;
        for (int t1 = 0; t1 < NPAIR; t1++)
            for (int t2 = t1 + 1; t2 < NPAIR; t2++)
                for (int c = 1; c < 256; c++) {
                    uint32_t masks[N];
                    for (int k = 0; k < N; k++)
                        masks[k] = MC[t1][1][k] ^ MC[t2][c][k];
                    int bad = count_bad(masks, 255, NULL);
                    if (bad < best) {
                        best = bad;
                        printf("new best bad=%d: x^%d + %d*x^%d\n", bad,
                               exponent(t1), c, exponent(t2));
                    }
                    if (bad == 0) hits++;
                }
        printf("binomials: best bad=%d, constant-rank-6 hits=%ld\n", best, hits);
        return 0;
    }
    if (!strcmp(mode, "trinomial")) {
        int best = 1 << 30;
        long hits = 0, total = 0;
        for (int t1 = 0; t1 < NPAIR; t1++)
            for (int t2 = t1 + 1; t2 < NPAIR; t2++)
                for (int t3 = t2 + 1; t3 < NPAIR; t3++)
                    for (int c2 = 1; c2 < 256; c2++)
                        for (int c3 = 1; c3 < 256; c3++) {
                            uint32_t masks[N];
                            for (int k = 0; k < N; k++)
                                masks[k] = MC[t1][1][k] ^ MC[t2][c2][k] ^
                                           MC[t3][c3][k];
                            int bad = count_bad(masks, maxbad, NULL);
                            total++;
                            if (bad <= maxbad) cnt[bad]++;
                            /* counts above maxbad are truncated: not exact */
                            if (bad <= maxbad && bad < best) {
                                best = bad;
                                printf("new best bad=%d: x^%d + %d*x^%d + %d*x^%d\n",
                                       bad, exponent(t1), c2, exponent(t2), c3,
                                       exponent(t3));
                                fflush(stdout);
                            }
                            if (bad == 0) hits++;
                        }
        if (best > maxbad)
            printf("trinomials: total=%ld, none with bad <= %d\n", total, maxbad);
        else
            printf("trinomials: total=%ld best bad=%d hits=%ld\n", total, best, hits);
        for (int b = 0; b <= maxbad && b < 256; b++)
            printf("  bad=%d: %ld\n", b, cnt[b]);
        return 0;
    }
    if (!strcmp(mode, "f2coef")) {
        int best = 1 << 30;
        long hits = 0;
        uint32_t bestA = 0;
        uint32_t masks[N] = {0}, A = 0;
        for (uint32_t g = 1; g < (1u << NPAIR); g++) {
            int t = __builtin_ctz(g);
            A ^= 1u << t;
            for (int k = 0; k < N; k++) masks[k] ^= MC[t][1][k];
            int bad = count_bad(masks, maxbad, NULL);
            if (bad <= maxbad) cnt[bad]++;
            if (bad <= maxbad && bad < best) { best = bad; bestA = A; }
            if (bad == 0) hits++;
        }
        printf("F2-coefficient DO: best bad=%d (A=0x%07x) hits=%ld\n", best,
               bestA, hits);
        for (int b = 0; b <= maxbad && b < 256; b++)
            printf("  bad=%d: %ld\n", b, cnt[b]);
        return 0;
    }
    fprintf(stderr, "unknown mode\n");
    return 1;
}
