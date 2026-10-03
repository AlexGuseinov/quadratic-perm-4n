#!/usr/bin/env python3
"""Task 1: rank-6 table (Pfaffian method) for all 2^28 alternating 8x8 forms over F_2,
computed bit-sliced with numpy via Pfaffians (no Gaussian elimination).

Mask m (28 bits): bit p <-> pair (i,j), i<j, lexicographic order.
rank 8  <=> Pf8(B) = 1
rank>=6 <=> some principal 6x6 Pfaffian = 1
rank 6  <=> Pf8 = 0 and OR_{28 six-subsets} Pf6 = 1

Output: rank6_pf.bin = 2^22 little-endian uint64 words, bit (m&63) of word (m>>6)
is 1 iff mask m has rank 6 (same layout as good8.bin).
"""
import os
DATA = os.environ.get('N8_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'))
import itertools, sys, time
import numpy as np

N = 8
PAIRS = [(i, j) for i in range(N) for j in range(i + 1, N)]
PIDX = {pr: p for p, pr in enumerate(PAIRS)}
assert len(PAIRS) == 28


def matchings(elems):
    """All perfect matchings of the (sorted) list elems, as lists of pairs."""
    if not elems:
        return [[]]
    first, rest = elems[0], elems[1:]
    out = []
    for k in range(len(rest)):
        partner = rest[k]
        remaining = rest[:k] + rest[k + 1:]
        for m in matchings(remaining):
            out.append([(first, partner)] + m)
    return out


M8 = matchings(list(range(8)))
assert len(M8) == 105
SUB6 = list(itertools.combinations(range(8), 6))
assert len(SUB6) == 28
M6 = {s: matchings(list(s)) for s in SUB6}
assert all(len(v) == 15 for v in M6.values())

LOWPAT = [np.uint64(0xAAAAAAAAAAAAAAAA), np.uint64(0xCCCCCCCCCCCCCCCC),
          np.uint64(0xF0F0F0F0F0F0F0F0), np.uint64(0xFF00FF00FF00FF00),
          np.uint64(0xFFFF0000FFFF0000), np.uint64(0xFFFFFFFF00000000)]


def rank6_chunk(w0, L):
    """Bitsliced rank-6 indicator for masks 64*w0 .. 64*(w0+L)-1 (L words)."""
    w = np.arange(w0, w0 + L, dtype=np.uint64)
    bits = []
    for p in range(28):
        if p < 6:
            bits.append(np.full(L, LOWPAT[p], dtype=np.uint64))
        else:
            b = (w >> np.uint64(p - 6)) & np.uint64(1)
            bits.append((np.uint64(0) - b).astype(np.uint64))  # 0 or all-ones
    x = {pr: bits[PIDX[pr]] for pr in PAIRS}
    # Pf8: XOR over the 105 matchings of AND of 4 pair bits
    pf8 = np.zeros(L, dtype=np.uint64)
    for m in M8:
        t = x[m[0]] & x[m[1]]
        t &= x[m[2]]
        t &= x[m[3]]
        pf8 ^= t
    # OR over 28 principal 6x6 Pfaffians
    any6 = np.zeros(L, dtype=np.uint64)
    for s in SUB6:
        pf6 = np.zeros(L, dtype=np.uint64)
        for m in M6[s]:
            t = x[m[0]] & x[m[1]]
            t &= x[m[2]]
            pf6 ^= t
        any6 |= pf6
    return any6 & ~pf8, pf8, any6


def popcount(a):
    return int(np.bitwise_count(a).sum())


def main():
    t0 = time.time()
    NW = 1 << 22
    CH = 1 << 18
    out = np.zeros(NW, dtype=np.uint64)
    n8 = 0
    nge6 = 0
    for w0 in range(0, NW, CH):
        r6, pf8, any6 = rank6_chunk(w0, CH)
        out[w0:w0 + CH] = r6
        n8 += popcount(pf8)
        nge6 += popcount(any6 | pf8)
    n6 = popcount(out)
    print(f"rank-8 forms: {n8}")
    print(f"rank>=6 forms: {nge6}  (rank-6 = {n6})")
    print(f"rank-6 forms: {n6}")
    # sanity: does Pf8=1 always imply some Pf6=1 (as it must: rank 8 => rank >= 6)?
    # (not needed for rank-6 definition, but a consistency check of the method)
    # theoretical count: q^{r(r-1)} prod_{i=0}^{2r-1}(q^{n-i}-1) / prod_{i=1}^r (q^{2i}-1)
    def cnt(n, r, q=2):
        num = q ** (r * (r - 1))
        for i in range(2 * r):
            num *= q ** (n - i) - 1
        den = 1
        for i in range(1, r + 1):
            den *= q ** (2 * i) - 1
        assert num % den == 0
        return num // den
    print(f"theoretical rank-6 count: {cnt(8,3)}, rank-8: {cnt(8,4)}")
    assert n6 == cnt(8, 3) and n8 == cnt(8, 4), "rank counts differ from the closed formula"
    out.astype('<u8').tofile(os.path.join(DATA, 'rank6_pf.bin'))
    print(f"time {time.time()-t0:.1f}s")
    # compare with good8.bin (the table of the first program), if it is present
    if not os.path.exists(os.path.join(DATA, 'good8.bin')):
        print("good8.bin not present: comparison with the table of the first program skipped")
        return
    g = np.fromfile(os.path.join(DATA, 'good8.bin'), dtype='<u8')
    print(f"good8.bin words: {len(g)} (expected {NW+1})")
    diff = g[:NW] ^ out
    nd = popcount(diff)
    print(f"differing bits vs good8.bin (first 2^22 words): {nd}")
    print("IDENTICAL" if nd == 0 else "DIFFERENT")
    if nd:
        idx = np.nonzero(diff)[0][:5]
        for wi in idx:
            for b in range(64):
                if (int(diff[wi]) >> b) & 1:
                    m = (int(wi) << 6) | b
                    print(f"  mask {m}: mine={(int(out[wi])>>b)&1} good8={(int(g[wi])>>b)&1}")
                    break
    print(f"extra word good8[2^22] = {int(g[NW]) if len(g) > NW else None}")


if __name__ == '__main__':
    main()
