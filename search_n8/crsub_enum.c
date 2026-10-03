/*
 * crsub_enum: exhaustive search for the component spaces of quadratic DU-4 permutations.
 *
 * Target: n-dimensional subspaces W of Alt(n, F_2) such that
 *   (i)  every nonzero element of W has rank n - 2, and
 *   (ii) for every point a != 0, P_a = {w in W : a in rad w} has dimension exactly 2.
 * By Charpin-Peng (L2 in NOTES_lemmas.md) the component space {B_b} of a quadratic DU-4
 * permutation of F_2^n satisfies (i) and (ii) (P_a = Im beta(a, .)^perp has dim s_a = 2).
 *
 * Search ("radical closure").  W_1 = <J>, J = e0^e1 + ... + e_{n-4}^e_{n-3}, radical
 * <e_{n-2}, e_{n-1}>; every solution is GL-equivalent to one containing J.
 * At a node W_d (d < n) a point r with dim(P_r cap W_d) <= 1 is chosen by a rule that
 * depends only on W_d.  In any solution W containing W_d, P_r has dimension 2, so W
 * contains an element c outside W_d with r in rad c; hence the coset c + W_d is "good"
 * (all its elements have rank n-2) and contains an element whose radical contains r.
 * The children are W_d + c for these cosets, one per orbit of a group H stabilizing W_d
 * and fixing r (H maps solutions containing W_d to solutions containing W_d).
 * So every solution is equivalent to one containing some leaf W_n; leaves are checked.
 * Pruning (all valid for solutions containing W_d):
 *   * dim(P_a cap W_d) <= 2 for every a;
 *   * a point a with dim(P_a cap W_d) = 1 needs >= 1 candidate coset, with dim 0 needs >= 3
 *     (the pencil P_a maps onto a 2-dim subspace of W / W_d);
 *   * at least 2^(n-d) - 1 good cosets.
 * Groups: G_1 = Stab(J) from explicit generators (checked); below that, point stabilizers
 * and coset stabilizers from random Schreier generators.  A proper subgroup only makes the
 * orbits finer (more work), never loses solutions.  Every group element used is checked at
 * run time to stabilize W_d (and to fix r where needed).
 *
 * Usage: crsub_enum n [--nosym] [--sample LEVEL K] [--seed S] [--out FILE] [--good FILE]
 *   --nosym        no symmetry reduction below level 1 (for validation)
 *   --sample L K   at levels >= L process at most K random children per node; node counts are
 *                  then weighted (Knuth estimator) to estimate the full tree
 *   --part I M     process only the level-3 nodes with index = I mod M (the random choices are
 *                  seeded by the path, so all parts see the same tree above level 3)
 *   --from K       skip level-3 nodes with index < K (resume)
 *   --list3        only list the level-3 nodes (index, hash and basis of W_3), for checking --part
 *   --nosym-from L no symmetry reduction at levels >= L (for comparison with other code)
 *   --cert FILE    write certificates for the orbit reductions at levels <= 2 (see expand)
 */
#include "crsub_common.h"
#include <time.h>
#include <math.h>

#define MAXGEN 48
#define MAXLV 9
#define NRAND 24          /* random-product stabilizer elements */
#define RANDLEV 3         /* use them at levels <= RANDLEV */

typedef struct { int ng; Mat g[MAXGEN]; Act a[MAXGEN]; } Group;

static Sub S[MAXLV + 1];
static uint32_t *Cl[MAXLV + 1];
static size_t nC[MAXLV + 1], capC[MAXLV + 1];
static uint64_t *Cb[MAXLV + 1];
static Group GRP[MAXLV + 1], PGRP[MAXLV + 1];
static uint8_t RT[256][4][128];
static int NOSYM_FROM = 99, NOSYM = 0, SAMPLE_L = 99, SAMPLE_K = 0;
static FILE *OUT, *CERT;
#define CERTLEV 2         /* certificates for the orbit reductions at levels <= CERTLEV */
static double tchild[MAXLV + 1], wnodes[MAXLV + 1], nodes[MAXLV + 1], tnode[MAXLV + 1], wcand[MAXLV + 1], wreps[MAXLV + 1];
static long long nsol = 0, prune_dim[MAXLV + 1], prune_cov[MAXLV + 1], prune_size[MAXLV + 1];
static uint32_t J;

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static inline int btest(const uint64_t *b, uint32_t i) { return (b[i >> 6] >> (i & 63)) & 1; }
static inline void bset(uint64_t *b, uint32_t i) { b[i >> 6] |= 1ull << (i & 63); }
static inline void bclr(uint64_t *b, uint32_t i) { b[i >> 6] &= ~(1ull << (i & 63)); }
static inline int inC(int d, uint32_t canon) { return canon && btest(Cb[d], sub_index(&S[d], canon)); }

static void push(int d, uint32_t x) {
    if (nC[d] == capC[d]) { capC[d] = capC[d] ? 2 * capC[d] : 1024; Cl[d] = realloc(Cl[d], capC[d] * 4); }
    Cl[d][nC[d]++] = x;
}

/* x(a, .) as an n-bit row, x a form */
static void build_rt(void) {
    for (int a = 1; a < (1 << n); a++) {
        uint8_t img[28];
        for (int p = 0; p < NP; p++) {
            int i = PI_[p], j = PJ_[p];
            img[p] = ((a >> i & 1) ? (1 << j) : 0) ^ ((a >> j & 1) ? (1 << i) : 0);
        }
        for (int c = 0; c < 4; c++) {
            RT[a][c][0] = 0;
            for (int v = 1; v < 128; v++) {
                int p = 7 * c + __builtin_ctz(v);
                RT[a][c][v] = RT[a][c][v & (v - 1)] ^ (p < NP ? img[p] : 0);
            }
        }
    }
}
static inline uint8_t rowf(int a, uint32_t x) {
    return RT[a][0][x & 127] ^ RT[a][1][(x >> 7) & 127] ^ RT[a][2][(x >> 14) & 127] ^ RT[a][3][(x >> 21) & 127];
}
static inline int matvec(Mat g, int v) {
    int w = 0;
    for (int i = 0; i < n; i++) w |= (__builtin_parity(g.r[i] & v)) << i;
    return w;
}

static void grp_clear(Group *G) { G->ng = 0; }
static int grp_add(Group *G, Mat g) {
    if (mat_is_id(g) || G->ng >= MAXGEN) return 0;
    for (int k = 0; k < G->ng; k++) if (mat_eq(G->g[k], g)) return 0;
    G->g[G->ng] = g; act_build(g, &G->a[G->ng]); G->ng++;
    return 1;
}
/* every generator maps W_d into W_d (and fixes r if r > 0) */
static void grp_check(const Group *G, int d, int r, const char *what) {
    for (int k = 0; k < G->ng; k++) {
        for (int t = 0; t < S[d].d; t++)
            if (sub_reduce(&S[d], act(&G->a[k], S[d].b[t]))) { fprintf(stderr, "BUG %s: generator does not stabilize W_%d\n", what, d); exit(2); }
        if (r > 0 && matvec(G->g[k], r) != r) { fprintf(stderr, "BUG %s: generator moves r\n", what); exit(2); }
    }
}

/* ---------- level 1 ---------- */
static void gens_level1(Group *G) {
    grp_clear(G);
    int m = n - 2;                               /* symplectic part on coords 0..m-1 */
    for (int v = 1; v < (1 << m); v++) {
        if (__builtin_popcount(v) > 2) continue;
        int phi = 0;                             /* J(v, x) = phi . x */
        for (int k = 0; k < m / 2; k++) {
            if (v >> (2 * k) & 1) phi |= 1 << (2 * k + 1);
            if (v >> (2 * k + 1) & 1) phi |= 1 << (2 * k);
        }
        Mat g = mat_id();
        for (int i = 0; i < n; i++) if (v >> i & 1) g.r[i] ^= phi;
        grp_add(G, g);
    }
    Mat g = mat_id(); g.r[n - 2] = 1 << (n - 1); g.r[n - 1] = 1 << (n - 2); grp_add(G, g);
    g = mat_id(); g.r[n - 2] ^= 1 << (n - 1); grp_add(G, g);
    for (int i = 0; i < m; i++) {
        g = mat_id(); g.r[n - 2] ^= 1 << i; grp_add(G, g);
        g = mat_id(); g.r[n - 1] ^= 1 << i; grp_add(G, g);
    }
    for (int k = 0; k < G->ng; k++) {
        Mat inv;
        if (!mat_inv(G->g[k], &inv) || act(&G->a[k], J) != J) { fprintf(stderr, "BUG: level-1 generator\n"); exit(2); }
    }
}

/* ---------- point data: for each a, a basis of {w(a, .) : w in W_d} (incremental) ---------- */
static uint8_t BAS[MAXLV + 1][256][8];
static int RK[MAXLV + 1][256];
static uint32_t *CAND[MAXLV + 1];
static size_t nCAND[MAXLV + 1], capCAND[MAXLV + 1];
static int FULLP[MAXLV + 1][256], nFULLP[MAXLV + 1];   /* points a with P_a inside W_d */

/* basis kept with distinct leading bits, sorted decreasing; v -> min(v, v ^ b) reduces */
static inline uint8_t reduce8(const uint8_t *b, int r, uint8_t v) {
    for (int s = 0; s < r; s++) { uint8_t t = v ^ b[s]; if (t < v) v = t; }
    return v;
}
static inline int in_span(int d, int a, uint8_t v) { return reduce8(BAS[d][a], RK[d][a], v) == 0; }
static void points_level1(void) {
    for (int a = 1; a < (1 << n); a++) {
        uint8_t v = rowf(a, J);
        RK[1][a] = 0;
        if (v) { BAS[1][a][0] = v; RK[1][a] = 1; }
    }
}
/* point data of W_d + c; returns 0 if some dim(P_a cap W_{d+1}) > 2 */
static int points_child(int d, uint32_t c) {
    int ok = 1, dd = S[d].d + 1;
    for (int a = 1; a < (1 << n); a++) {
        int r = RK[d][a];
        memcpy(BAS[d + 1][a], BAS[d][a], 8);
        uint8_t v = reduce8(BAS[d][a], r, rowf(a, c));
        if (v) {
            uint8_t *b = BAS[d + 1][a];
            b[r] = v;
            for (int s = r; s > 0 && b[s] > b[s - 1]; s--) { uint8_t t = b[s]; b[s] = b[s - 1]; b[s - 1] = t; }
            r++;
        }
        RK[d + 1][a] = r;
        if (dd - r > 2) ok = 0;
    }
    nFULLP[d + 1] = 0;
    for (int a = 1; a < (1 << n); a++) if (dd - RK[d + 1][a] == 2) FULLP[d + 1][nFULLP[d + 1]++] = a;
    return ok;
}

static void cand_push(int d, uint32_t x) {
    if (nCAND[d] == capCAND[d]) { capCAND[d] = capCAND[d] ? 2 * capCAND[d] : 1024; CAND[d] = realloc(CAND[d], capCAND[d] * 4); }
    CAND[d][nCAND[d]++] = x;
}

/* returns r > 0, or 0 if the node is pruned, or -1 if d == n and the leaf is a solution */
static int choose_point(int d) {
    int N = 1 << n;
    if (d == n) {
        for (int a = 1; a < N; a++) if (d - RK[d][a] != 2) { prune_dim[d]++; return 0; }
        return -1;
    }
    if (nC[d] < (size_t)((1 << (n - d)) - 1)) { prune_size[d]++; return 0; }
    int r = 0;
    if (d <= 2) {
        for (int a = 1; a < N && !r; a++) if (d - RK[d][a] == 1) r = a;
        for (int a = 1; a < N && !r; a++) if (d - RK[d][a] == 0) r = a;
    } else {
        long best = -1;
        for (int a = 1; a < N; a++) {
            int dp = d - RK[d][a];
            if (dp == 2) continue;
            long cnt = 0;
            for (size_t t = 0; t < nC[d]; t++) cnt += in_span(d, a, rowf(a, Cl[d][t]));
            if (cnt == 0 || (dp == 0 && cnt < 3)) { prune_cov[d]++; return 0; }
            long key = cnt * 4 + (dp == 1 ? 0 : 1);
            if (best < 0 || key < best) { best = key; r = a; }
        }
    }
    nCAND[d] = 0;
    for (size_t t = 0; t < nC[d]; t++) if (in_span(d, r, rowf(r, Cl[d][t]))) cand_push(d, Cl[d][t]);
    if (nCAND[d] == 0) { prune_cov[d]++; return 0; }
    return r;
}

/* ---------- child sets ---------- */
/* small hash sets (fast membership) for the sets of levels >= HLEV */
#define HLEV 4
static uint32_t *HT[MAXLV + 1]; static size_t HM[MAXLV + 1], HCAP[MAXLV + 1];
static void hs_build(int d) {
    size_t sz = 64; while (sz < 4 * nC[d]) sz <<= 1;
    if (sz > HCAP[d]) { free(HT[d]); HT[d] = calloc(sz, 4); HCAP[d] = sz; }
    HM[d] = sz - 1;
    for (size_t t = 0; t < nC[d]; t++) {
        uint32_t k = Cl[d][t]; size_t i = (k * 2654435761u) & HM[d];
        while (HT[d][i]) i = (i + 1) & HM[d];
        HT[d][i] = k;
    }
}
static inline int hs_has(int d, uint32_t k) {
    size_t i = (k * 2654435761u) & HM[d];
    while (HT[d][i]) { if (HT[d][i] == k) return 1; i = (i + 1) & HM[d]; }
    return 0;
}
static void hs_clear(int d) { memset(HT[d], 0, (HM[d] + 1) * 4); }

static void make_child(int d, uint32_t c) {
    S[d + 1] = S[d];
    sub_add(&S[d + 1], c);
    nC[d + 1] = 0;
    for (size_t t = 0; t < nC[d]; t++) {
        uint32_t y = Cl[d][t];
        uint32_t z = sub_reduce(&S[d], y ^ c);
        if (!z) continue;
        if (d >= HLEV ? !hs_has(d, z) : !btest(Cb[d], sub_index(&S[d], z))) continue;
        uint32_t y2 = sub_reduce(&S[d + 1], y);
        uint32_t i2 = sub_index(&S[d + 1], y2);
        if (btest(Cb[d + 1], i2)) continue;
        /* filter: if P_a is inside W_{d+1}, no element outside W_{d+1} has a in its radical */
        int bad = 0;
        for (int k = 0; k < nFULLP[d + 1] && !bad; k++) {
            int a = FULLP[d + 1][k];
            if (in_span(d + 1, a, rowf(a, y2))) bad = 1;
        }
        if (bad) continue;
        bset(Cb[d + 1], i2);
        push(d + 1, y2);
    }
    if (d + 1 >= HLEV) hs_build(d + 1);
}
static void clear_child(int d) {
    for (size_t t = 0; t < nC[d]; t++) bclr(Cb[d], sub_index(&S[d], Cl[d][t]));
    if (d >= HLEV) hs_clear(d);
    nC[d] = 0;
}

/* ---------- random group elements (random words of length 64) ---------- */
static Mat rand_elem(const Group *H) {
    Mat g = mat_id();
    for (int t = 0; t < 64; t++) g = mat_mul(g, H->g[rnd() % H->ng]);
    return g;
}

/* ---------- stabilizer of a point ---------- */
static void stab_point(const Group *H, int r, Group *out) {
    grp_clear(out);
    if (H->ng == 0) return;
    static Mat u[256]; static int seen[256], orb[256];
    memset(seen, 0, sizeof seen);
    int len = 0; orb[len++] = r; seen[r] = 1; u[r] = mat_id();
    for (int h = 0; h < len; h++) {
        int v = orb[h];
        for (int k = 0; k < H->ng; k++) {
            int w = matvec(H->g[k], v);
            if (!seen[w]) { seen[w] = 1; orb[len++] = w; u[w] = mat_mul(H->g[k], u[v]); }
        }
    }
    /* nearly uniform elements g of H, corrected by the transversal: u_{g r}^-1 g fixes r */
    for (int att = 0; att < 2 * NRAND && out->ng < NRAND; att++) {
        Mat g = rand_elem(H), inv;
        mat_inv(u[matvec(g, r)], &inv);
        grp_add(out, mat_mul(inv, g));
    }
    for (int att = 0; att < 6 * MAXGEN && out->ng < MAXGEN - 8; att++) {
        int v = orb[rnd() % len], k = rnd() % H->ng;
        int w = matvec(H->g[k], v);
        Mat inv; mat_inv(u[w], &inv);
        grp_add(out, mat_mul(inv, mat_mul(H->g[k], u[v])));   /* u_w^-1 s u_v fixes r */
    }
}

/* ---------- orbits on candidate cosets ---------- */
typedef struct { uint32_t *key; int *val; size_t mask; } Hash;
static void hash_init(Hash *h, size_t m) {
    size_t sz = 16; while (sz < 2 * m) sz <<= 1;
    h->key = malloc(sz * 4); h->val = malloc(sz * sizeof(int)); h->mask = sz - 1;
    memset(h->key, 0, sz * 4);
}
static void hash_free(Hash *h) { free(h->key); free(h->val); }
static inline void hash_put(Hash *h, uint32_t k, int v) {
    size_t i = (k * 2654435761u) & h->mask;
    while (h->key[i]) i = (i + 1) & h->mask;
    h->key[i] = k; h->val[i] = v;
}
static inline int hash_get(const Hash *h, uint32_t k) {
    size_t i = (k * 2654435761u) & h->mask;
    while (h->key[i]) { if (h->key[i] == k) return h->val[i]; i = (i + 1) & h->mask; }
    return -1;
}

static void dfs(int d, double w, uint64_t seed);
static uint64_t w_hash(int d);
static uint64_t mix64(uint64_t z) { z += 0x9E3779B97F4A7C15ull; z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }

static void expand(int d, int r, double w, double t0, uint64_t seed) {
    const Group *H = &PGRP[d];
    size_t m = nCAND[d];
    uint32_t *cand = CAND[d];
    /* orbit reps (with BFS data for stabilizers) */
    int *order = NULL, *ostart = NULL, *olen = NULL, nor = 0;
    Mat *u = NULL; Hash hsh; int have_hash = 0;
    if (H->ng == 0) {
        nor = (int)m;
    } else {
        hash_init(&hsh, m); have_hash = 1;
        for (size_t t = 0; t < m; t++) hash_put(&hsh, cand[t], (int)t);
        order = malloc(m * sizeof(int)); ostart = malloc((m + 1) * sizeof(int)); olen = malloc((m + 1) * sizeof(int));
        u = malloc(m * sizeof(Mat));
        char *vis = calloc(m, 1);
        int len = 0;
        for (size_t t = 0; t < m; t++) {
            if (vis[t]) continue;
            ostart[nor] = len; vis[t] = 1; order[len++] = (int)t; u[t] = mat_id();
            for (int h = ostart[nor]; h < len; h++) {
                int x = order[h];
                for (int k = 0; k < H->ng; k++) {
                    uint32_t y = sub_reduce(&S[d], act(&H->a[k], cand[x]));
                    int iy = hash_get(&hsh, y);
                    if (iy < 0) { fprintf(stderr, "BUG: group does not preserve candidates (d=%d)\n", d); exit(2); }
                    if (!vis[iy]) { vis[iy] = 1; order[len++] = iy; u[iy] = mat_mul(u[x], H->g[k]); }
                }
            }
            olen[nor] = len - ostart[nor];
            nor++;
        }
        free(vis);
    }
    wcand[d] += w * m;
    if (CERT && d <= CERTLEV) {
        /* NODE line: level, basis of W_d, point r, orbit representatives; then one line per
           candidate x: x, index of its orbit representative, matrix u (rows, hex) with
           rep^u = x (mod W_d), where B^u(x, y) = B(u x, u y) */
        fprintf(CERT, "NODE %d %016llx", d, (unsigned long long)w_hash(d));
        for (int k = 0; k < S[d].d; k++) fprintf(CERT, " %u", S[d].b[k]);
        fprintf(CERT, " r %d reps %d", r, nor);
        for (int o = 0; o < nor; o++) fprintf(CERT, " %u", H->ng == 0 ? cand[o] : cand[order[ostart[o]]]);
        fprintf(CERT, "\n");
        if (H->ng > 0)
            for (int o = 0; o < nor; o++)
                for (int h = ostart[o]; h < ostart[o] + olen[o]; h++) {
                    int x = order[h];
                    fprintf(CERT, "%u %d ", cand[x], o);
                    for (int i = 0; i < n; i++) fprintf(CERT, "%02x", u[x].r[i]);
                    fprintf(CERT, "\n");
                }
    }
    wreps[d] += w * nor;
    /* which children to process */
    int *pick = malloc((nor + 1) * sizeof(int)), npick = nor;
    for (int i = 0; i < nor; i++) pick[i] = i;
    double cw = w;
    if (d >= SAMPLE_L && SAMPLE_K > 0 && nor > SAMPLE_K) {
        for (int i = 0; i < SAMPLE_K; i++) { int j = i + rnd() % (nor - i); int tt = pick[i]; pick[i] = pick[j]; pick[j] = tt; }
        npick = SAMPLE_K; cw = w * (double)nor / SAMPLE_K;
    }
    tnode[d] += now() - t0;
    for (int pi = 0; pi < npick; pi++) {
        int o = pick[pi];
        rng_s = mix64(seed ^ ((uint64_t)o * 0x9E3779B97F4A7C15ull + 7)) | 1;   /* path-determined */
        uint32_t c = (H->ng == 0) ? cand[o] : cand[order[ostart[o]]];
        /* child group: stabilizer of the coset c + W_d in H */
        Group *CG = &GRP[d + 1];
        grp_clear(CG);
        if (H->ng > 0 && olen[o] >= 1) {
            int st = ostart[o], L = olen[o];
            if (d <= RANDLEV)
                for (int att = 0; att < 2 * NRAND && CG->ng < NRAND; att++) {
                    Mat g = rand_elem(H), inv;
                    int ix = hash_get(&hsh, sub_reduce(&S[d], act_slow(g, c)));
                    if (ix < 0) { fprintf(stderr, "BUG: random element leaves the candidates\n"); exit(2); }
                    mat_inv(u[ix], &inv);
                    Mat h = mat_mul(g, inv);                   /* rep^(g u_x^-1) = rep */
                    if (sub_reduce(&S[d], act_slow(h, c)) != c) { fprintf(stderr, "BUG: random stabilizer element\n"); exit(2); }
                    grp_add(CG, h);
                }
            for (int att = 0; att < 6 * MAXGEN && CG->ng < MAXGEN - 8; att++) {
                int x = order[st + rnd() % L], k = rnd() % H->ng;
                uint32_t y = sub_reduce(&S[d], act(&H->a[k], cand[x]));
                int iy = hash_get(&hsh, y);
                Mat inv; mat_inv(u[iy], &inv);
                Mat g = mat_mul(mat_mul(u[x], H->g[k]), inv);   /* rep^g = rep */
                if (sub_reduce(&S[d], act_slow(g, c)) != c) { fprintf(stderr, "BUG: Schreier generator\n"); exit(2); }
                grp_add(CG, g);
            }
        }
        double tc = now();
        if (!points_child(d, c)) {                   /* child pruned before building its set */
            nodes[d + 1] += 1; wnodes[d + 1] += cw; prune_dim[d + 1]++;
            tnode[d + 1] += now() - tc;
            continue;
        }
        make_child(d, c);
        tchild[d + 1] += now() - tc;
        dfs(d + 1, cw, mix64(seed ^ ((uint64_t)o * 0x2545F4914F6CDD1Dull + 1)));
        clear_child(d + 1);
    }
    free(pick);
    if (have_hash) { hash_free(&hsh); free(order); free(ostart); free(olen); free(u); }
}

static int PART_I = 0, PART_M = 1, LIST3 = 0; static long PART_FROM = 0, L3IDX = -1;
static uint64_t w_hash(int d) {                   /* hash of the reduced echelon basis of W_d */
    uint32_t b[32]; int k = S[d].d;
    memcpy(b, S[d].b, 4 * k);
    for (int i = 0; i < k; i++) for (int j = i + 1; j < k; j++) if (b[j] > b[i]) { uint32_t t = b[i]; b[i] = b[j]; b[j] = t; }
    uint64_t h = 1469598103934665603ull;
    for (int i = 0; i < k; i++) h = mix64(h ^ b[i]);
    return h;
}
#define PARTLEV 3
static double T_START;
static void dfs(int d, double w, uint64_t seed) {
    if (d == PARTLEV) {
        L3IDX++;
        if (L3IDX % PART_M != PART_I || L3IDX < PART_FROM) { tchild[d] = 0; return; }
        if (LIST3) {
            fprintf(stderr, "L3 %ld W3 %016llx basis", L3IDX, (unsigned long long)w_hash(d));
            for (int k = 0; k < S[d].d; k++) fprintf(stderr, " %u", S[d].b[k]);
            fprintf(stderr, "\n"); tchild[d] = 0; return;
        }
    }
    rng_s = seed ? seed : 1;
    double t0 = now() - tchild[d];
    tchild[d] = 0;
    nodes[d] += 1; wnodes[d] += w;
    grp_check(&GRP[d], d, 0, "node");
    int r = choose_point(d);
    if (r == -1) {
        nsol++;
        fprintf(OUT, "SOLUTION");
        for (int k = 0; k < S[d].d; k++) fprintf(OUT, " %u", S[d].b[k]);
        fprintf(OUT, "\n"); fflush(OUT);
        tnode[d] += now() - t0;
        return;
    }
    if (r == 0) {
        tnode[d] += now() - t0;
        if (d == PARTLEV && n >= 8) {
            fprintf(stderr, "L3 %ld W3 %016llx pruned: solutions %lld, %.0fs\n", L3IDX, (unsigned long long)w_hash(d),
                    nsol, now() - T_START);
            fflush(stderr);
        }
        return;
    }
    if (NOSYM || d >= NOSYM_FROM || GRP[d].ng == 0) grp_clear(&PGRP[d]);
    else stab_point(&GRP[d], r, &PGRP[d]);
    grp_check(&PGRP[d], d, r, "point stabilizer");
    expand(d, r, w, t0, seed);
    if (d == PARTLEV && n >= 8) {
        fprintf(stderr, "L3 %ld W3 %016llx done: solutions %lld, nodes4 %.0f nodes5 %.0f nodes6 %.0f nodes7 %.0f, %.0fs\n",
                L3IDX, (unsigned long long)w_hash(d), nsol, nodes[4], nodes[5], nodes[6], nodes[7], now() - T_START);
        fflush(stderr);
    }
}

/* self-test of the table-driven primitives against direct computations */
static void self_test(void) {
    for (int t = 0; t < 20000; t++) {
        uint32_t x = (uint32_t)rnd() & FULL;
        if (good(x) != (rank_mask(x) == n - 2)) { fprintf(stderr, "SELFTEST: GOOD\n"); exit(2); }
        Mat g; do { for (int i = 0; i < n; i++) g.r[i] = rnd() & ((1 << n) - 1); } while (!mat_inv(g, &(Mat){0}));
        Act A; act_build(g, &A);
        uint32_t y = act(&A, x);
        /* direct: y(e_k, e_l) = x(g e_k, g e_l) */
        uint32_t z = 0;
        for (int p = 0; p < NP; p++) {
            int k = PI_[p], l = PJ_[p], gk = 0, gl = 0;
            for (int i = 0; i < n; i++) { gk |= (g.r[i] >> k & 1) << i; gl |= (g.r[i] >> l & 1) << i; }
            int v = 0;
            for (int q = 0; q < NP; q++) if (x >> q & 1) {
                int i = PI_[q], j = PJ_[q];
                v ^= ((gk >> i & 1) & (gl >> j & 1)) ^ ((gk >> j & 1) & (gl >> i & 1));
            }
            z |= (uint32_t)v << p;
        }
        if (y != z || act_slow(g, x) != z) { fprintf(stderr, "SELFTEST: action\n"); exit(2); }
        int a = 1 + (int)(rnd() % ((1 << n) - 1)), row = 0;
        for (int l = 0; l < n; l++) {
            int v = 0;
            for (int q = 0; q < NP; q++) if (x >> q & 1) {
                int i = PI_[q], j = PJ_[q];
                v ^= ((a >> i & 1) & (j == l)) ^ ((a >> j & 1) & (i == l));
            }
            row |= v << l;
        }
        if (rowf(a, x) != row) { fprintf(stderr, "SELFTEST: rowf\n"); exit(2); }
    }
    fprintf(stderr, "self-test passed\n");
}

int main(int argc, char **argv) {
    init_pairs(atoi(argv[1]));
    const char *outp = NULL, *goodp = NULL;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "--nosym")) NOSYM = 1;
        else if (!strcmp(argv[i], "--sample")) { SAMPLE_L = atoi(argv[++i]); SAMPLE_K = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--seed")) rng_s ^= strtoull(argv[++i], 0, 10) * 0x9E3779B97F4A7C15ull;
        else if (!strcmp(argv[i], "--out")) outp = argv[++i];
        else if (!strcmp(argv[i], "--good")) goodp = argv[++i];
        else if (!strcmp(argv[i], "--part")) { PART_I = atoi(argv[++i]); PART_M = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--from")) PART_FROM = atol(argv[++i]);
        else if (!strcmp(argv[i], "--list3")) LIST3 = 1;
        else if (!strcmp(argv[i], "--nosym-from")) NOSYM_FROM = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--cert")) CERT = fopen(argv[++i], "w");
    }
    OUT = outp ? fopen(outp, "w") : stdout;
    long long hist[9];
    double T0 = now();
    build_good_cached(goodp, hist);
    build_rt();
    { uint64_t keep = rng_s; self_test(); rng_s = keep; }   /* the search does not depend on the test */
    J = 0;
    for (int k = 0; k < (n - 2) / 2; k++) J |= 1u << PAIR[2 * k][2 * k + 1];
    fprintf(stderr, "n=%d NP=%d GOOD ready (%.1fs)\n", n, NP, now() - T0);
    for (int d = 1; d <= n; d++) Cb[d] = calloc(((size_t)1 << (NP - d)) / 64 + 1, 8);
    sub_init(&S[1]); sub_add(&S[1], J);
    for (uint32_t i = 1; i < (1u << (NP - 1)); i++) {
        uint32_t rep = sub_unindex(&S[1], i);
        if (good(rep) && good(rep ^ J)) { bset(Cb[1], i); push(1, rep); }
    }
    fprintf(stderr, "level 1: %zu good cosets (%.1fs)\n", nC[1], now() - T0);
    points_level1();
    gens_level1(&GRP[1]);
    if (NOSYM) grp_clear(&GRP[1]);
    T_START = now();
    dfs(1, 1.0, mix64(rng_s));
    fprintf(stderr, "done in %.1fs, solutions (leaves passing (ii)): %lld\n", now() - T0, nsol);
    fprintf(stderr, "level  nodes      est.nodes     est.cands     est.reps      s/node    prune(dim,cov,size)\n");
    double tot = 0;
    for (int d = 1; d <= n; d++) {
        double per = nodes[d] > 0 ? tnode[d] / nodes[d] : 0;
        tot += per * wnodes[d];
        fprintf(stderr, "%3d %10.0f %13.4g %13.4g %13.4g %10.3g   %lld %lld %lld\n", d, nodes[d], wnodes[d], wcand[d], wreps[d],
                per, prune_dim[d], prune_cov[d], prune_size[d]);
    }
    fprintf(stderr, "COUNTS");
    for (int d = 1; d <= n; d++) fprintf(stderr, " n%d=%.0f", d, nodes[d]);
    for (int d = 1; d <= n; d++) fprintf(stderr, " pd%d=%lld pc%d=%lld ps%d=%lld", d, prune_dim[d], d, prune_cov[d], d, prune_size[d]);
    fprintf(stderr, "\n");
    fprintf(stderr, "estimated total node time: %.4g s (%.4g CPU-hours)\n", tot, tot / 3600);
    if (OUT != stdout) fclose(OUT);
    if (CERT) fclose(CERT);
    return 0;
}
