"""Shared primitives for the second implementation (Task 2).

Conventions (from the specification):
  point a: integer 1..255, bit i = coordinate i of F_2^8
  form mask m: 28 bits, bit p <-> pair (i,j) (i<j, lexicographic):
      B(x,y) = sum_{p=(i,j) in m} (x_i y_j + x_j y_i) mod 2
  rank-6 table: bitset, bit (m&63) of word (m>>6)   (produced by task1_pfaffian.py)

Key linear-algebra fact used throughout (derived from the definition of B):
  a in rad(m)  <=>  phi_a(m) := M_m a = 0, where (M_m a)_i = sum_j B_ij a_j.
  phi_a is LINEAR in the mask m.  For pair p=(i,j): phi_a(e_p) = a_j*e_i + a_i*e_j.
"""
import os
DATA = os.environ.get('N8_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'))
import numpy as np

NPAIR = 28
PAIRS = [(i, j) for i in range(8) for j in range(i + 1, 8)]
assert len(PAIRS) == NPAIR

RANK6_FILE = os.path.join(DATA, 'rank6_pf.bin')
_GOOD = None


def good_table():
    global _GOOD
    if _GOOD is None:
        _GOOD = np.fromfile(RANK6_FILE, dtype='<u8')
        assert len(_GOOD) == 1 << 22
    return _GOOD


def is_good(m):
    """m: integer array of masks -> bool array (rank exactly 6), from the Pfaffian table."""
    G = good_table()
    m = np.asarray(m)
    w = G[(m >> 6).astype(np.int64)]
    sh = (m & 63).astype(np.uint64)
    return ((w >> sh) & np.uint64(1)).astype(bool)


# ---------------------------------------------------------------------------
# phi tables:  PHI_T[k][b] is a length-256 uint8 row over points a:
#   PHI_T[k][b][a] = phi_a(b << 8k)
# so that phi_all(Y)[n, a] = XOR_k PHI_T[k][(Y[n] >> 8k) & 255][a]
# ---------------------------------------------------------------------------
def _build_phi():
    unit = np.zeros((NPAIR, 256), dtype=np.uint8)  # unit[p][a] = phi_a(e_p)
    for p, (i, j) in enumerate(PAIRS):
        for a in range(256):
            unit[p, a] = (((a >> j) & 1) << i) | (((a >> i) & 1) << j)
    PHI_T = np.zeros((4, 256, 256), dtype=np.uint8)
    for k in range(4):
        for b in range(256):
            row = np.zeros(256, dtype=np.uint8)
            for t in range(8):
                p = 8 * k + t
                if (b >> t) & 1 and p < NPAIR:
                    row ^= unit[p]
            PHI_T[k, b] = row
    return unit, PHI_T


UNIT, PHI_T = _build_phi()


def phi_all(Y):
    """Y: 1-D int array of masks -> (len(Y), 256) uint8, [n, a] = phi_a(Y[n])."""
    Y = np.asarray(Y).astype(np.int64)
    out = PHI_T[0][Y & 255]
    out ^= PHI_T[1][(Y >> 8) & 255]
    out ^= PHI_T[2][(Y >> 16) & 255]
    out ^= PHI_T[3][(Y >> 24) & 255]
    return out


def bilinear(m, x, y):
    """Direct evaluation of B(x,y) from the definition (scalar, for tests)."""
    s = 0
    for p, (i, j) in enumerate(PAIRS):
        if (m >> p) & 1:
            s ^= (((x >> i) & 1) & ((y >> j) & 1)) ^ (((x >> j) & 1) & ((y >> i) & 1))
    return s


def radical_points_direct(m):
    """Scalar: radical points of m by the definition B(a,y)=0 for all y (tests only)."""
    return [a for a in range(1, 256) if all(bilinear(m, a, 1 << t) == 0 for t in range(8))]


# ---------------------------------------------------------------------------
# bit-set helpers for 256-bit point sets stored as 4 uint64 words (little endian)
# ---------------------------------------------------------------------------
POINTMASK = np.array([~np.uint64(1), ~np.uint64(0), ~np.uint64(0), ~np.uint64(0)], dtype=np.uint64)


def pack_bool256(U):
    """(n,256) bool -> (n,4) uint64 bitsets (bit a <-> point a); point 0 cleared."""
    P = np.packbits(U, axis=1, bitorder='little')  # (n,32) uint8
    P = np.ascontiguousarray(P).view('<u8').reshape(-1, 4).copy()
    P &= POINTMASK
    return P


def unpack256(I):
    """(n,4) uint64 -> (n,256) uint8 0/1."""
    I = np.ascontiguousarray(I, dtype='<u8')
    return np.unpackbits(I.view(np.uint8).reshape(-1, 32), axis=1, bitorder='little')


def points_to_mask(pts_bool256):
    """(256,) bool -> (4,) uint64 bitset."""
    return pack_bool256(np.asarray(pts_bool256, dtype=bool).reshape(1, 256))[0]


# ---------------------------------------------------------------------------
# linear algebra on forms (28-bit masks)
# ---------------------------------------------------------------------------
def span(basis):
    els = [0]
    for b in basis:
        els = els + [e ^ b for e in els]
    return els


def rank_gf2(vecs):
    rows = [v for v in vecs]
    r = 0
    piv_rows = []
    for v in rows:
        for pv in piv_rows:
            v = min(v, v ^ pv)
        if v:
            piv_rows.append(v)
            r += 1
    return r


def echelon(basis):
    """Reduced echelon basis: returns (vecs, pivots) with pivot = highest bit of each
    vector, and every vector having 0 at the other vectors' pivots."""
    vecs = []
    for b in basis:
        v = b
        for w in vecs:
            if (v >> (w.bit_length() - 1)) & 1:
                v ^= w
        assert v != 0, "basis not independent"
        p = v.bit_length() - 1
        vecs = [w ^ v if (w >> p) & 1 else w for w in vecs]
        vecs.append(v)
    pivs = [v.bit_length() - 1 for v in vecs]
    for v in vecs:
        for q in pivs:
            if q != v.bit_length() - 1:
                assert not (v >> q) & 1
    return vecs, pivs


def compress(Y, pivs_desc):
    """Remove the (zero) bits at positions pivs_desc (sorted descending) from the masks Y."""
    for p in pivs_desc:
        low = (1 << p) - 1
        Y = (Y & low) | ((Y >> (p + 1)) << p)
    return Y


def expand(Z, pivs_asc):
    """Inverse of compress: insert zero bits at positions pivs_asc (sorted ascending)."""
    for p in pivs_asc:
        low = (1 << p) - 1
        Z = (Z & low) | ((Z >> p) << (p + 1))
    return Z


def dimP_direct(basis):
    """dimP_W(a) for all a (array of 256; index 0 unused) from radicals of all elements."""
    els = np.array([e for e in span(basis) if e != 0], dtype=np.int64)
    Z = (phi_all(els) == 0)  # (2^d-1, 256)
    cnt = Z.sum(0).astype(np.int64)
    cnt[0] = 0
    dimP = np.zeros(256, dtype=np.int64)
    for a in range(1, 256):
        c = int(cnt[a]) + 1
        assert c & (c - 1) == 0, "P_a(W) not a subspace?!"
        dimP[a] = c.bit_length() - 1
    return dimP


def coset_incidence_direct(Y, Wels):
    """For coset reps Y (array) and all elements Wels of W: (n,4) bitsets of the union of
    radicals over the coset elements (point a in set <=> coset meets R_a). Direct: phi_a == 0."""
    n = len(Y)
    out = np.zeros((n, 4), dtype=np.uint64)
    CH = 1 << 15
    Y = np.asarray(Y).astype(np.int64)
    for s in range(0, n, CH):
        y = Y[s:s + CH]
        U = np.zeros((len(y), 256), dtype=bool)
        for w in Wels:
            U |= (phi_all(y ^ w) == 0)
        out[s:s + CH] = pack_bool256(U)
    return out


def brute_C(basis, pivs):
    """C(W) from scratch by the definition, over ALL 2^(28-d) cosets.
    Canonical reps: zero on the pivot positions `pivs` (must be a valid pivot set for W).
    Returns (sorted reps, incidence bitsets, dimP)."""
    d = len(basis)
    vecs, pv = echelon(basis)
    Wels = span(basis)
    dimP = dimP_direct(basis)
    full = np.zeros(256, dtype=bool)
    full[1:] = dimP[1:] == 2
    fullmask = points_to_mask(full)
    pa = sorted(pivs)
    # check that pivs is a valid pivot set: projection of W onto pivs coords is bijective
    proj = set()
    for e in Wels:
        proj.add(tuple((e >> p) & 1 for p in pa))
    assert len(proj) == len(Wels), "invalid pivot set"
    nrep = 1 << (28 - d)
    reps_all = []
    CH = 1 << 22
    for s in range(0, nrep, CH):
        z = np.arange(s, min(nrep, s + CH), dtype=np.int64)
        y = expand(z, pa)
        ok = np.ones(len(y), dtype=bool)
        for w in Wels:
            ok &= is_good(y ^ w)
        reps_all.append(y[ok])
    Y = np.concatenate(reps_all)
    I = coset_incidence_direct(Y, Wels)
    keep = ((I & fullmask) == 0).all(1)
    Y, I = Y[keep], I[keep]
    order = np.argsort(Y)
    return Y[order], I[order], dimP
