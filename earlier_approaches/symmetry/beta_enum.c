/* Stage 1 of the symmetric search (n <= 8): enumerate every beta in a given
 * F_2-linear space (spanned by basis vectors read from stdin) and keep those with
 *   rank beta(a, .) = n-2 for every a != 0   (s_a = 2)   and
 *   rank B_b = n-2 for every b != 0          (all components of rank n-2).
 * Input:  n dim [mode], then one offset line and dim basis lines, each of n(n-1)/2 hex
 *         bytes c_ij = beta(e_i, e_j) for i < j in lexicographic pair order
 *         (0,1),(0,2),...,(n-2,n-1); the space enumerated is offset + span(basis).
 *         An optional third number on the first line, mode: 0 (default) = both
 *         conditions; 1 = only the component ranks (constant-rank test, Conjecture H).
 * Output: one line per surviving beta (hex bytes), then "count N of 2^dim".
 * Build:  cc -O3 -o beta_enum beta_enum.c
 */
#include <stdio.h>
#include <stdint.h>

static int rank8(uint8_t *r) {
    uint8_t m[8];
    int rk = 0;
    for (int i = 0; i < 8; i++) m[i] = r[i];
    for (int bit = 7; bit >= 0; bit--) {
        int p = -1;
        for (int i = rk; i < 8; i++) if (m[i] >> bit & 1) { p = i; break; }
        if (p < 0) continue;
        uint8_t t = m[p]; m[p] = m[rk]; m[rk] = t;
        for (int i = 0; i < 8; i++) if (i != rk && (m[i] >> bit & 1)) m[i] ^= m[rk];
        rk++;
    }
    return rk;
}

int main(void) {
    int n, dim, mode = 0;
    if (scanf("%d %d", &n, &dim) != 2 || dim > 40 || n > 8 || n < 2) return 1;
    {   /* optional mode on the same line */
        int ch;
        while ((ch = getchar()) == ' ') ;
        if (ch >= '0' && ch <= '9') mode = ch - '0'; else if (ch != EOF) ungetc(ch, stdin);
    }
    const int N = 1 << n;
    static uint8_t basis[40][8][8];
    uint8_t c[8][8] = {{0}};
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            unsigned v;
            if (scanf("%2x", &v) != 1) return 1;
            c[i][j] = c[j][i] = (uint8_t)v;
        }
    for (int k = 0; k < dim; k++) {
        for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) basis[k][i][j] = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) {
                unsigned v;
                if (scanf("%2x", &v) != 1) return 1;
                basis[k][i][j] = basis[k][j][i] = (uint8_t)v;
            }
    }
    uint64_t total = 1ULL << dim, count = 0;
    for (uint64_t t = 0; t < total; t++) {
        if (t) {
            int k = __builtin_ctzll(t);
            for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) c[i][j] ^= basis[k][i][j];
        }
        int ok = 1;
        for (int a = 1; a < N && ok && mode == 0; a++) {
            uint8_t rows[8] = {0};
            for (int j = 0; j < n; j++) {
                uint8_t v = 0;
                for (int i = 0; i < n; i++) if (a >> i & 1) v ^= c[i][j];
                rows[j] = v;
            }
            if (rank8(rows) != n - 2) ok = 0;
        }
        for (int b = 1; b < N && ok; b++) {
            uint8_t rows[8] = {0};
            for (int i = 0; i < n; i++) {
                uint8_t v = 0;
                for (int j = 0; j < n; j++) v |= (uint8_t)((__builtin_popcount(b & c[i][j]) & 1) << j);
                rows[i] = v;
            }
            if (rank8(rows) != n - 2) ok = 0;
        }
        if (ok) {
            count++;
            for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) printf("%02x", c[i][j]);
            printf("\n");
        }
    }
    printf("count %llu of 2^%d\n", (unsigned long long)count, dim);
    return 0;
}
