#!/usr/bin/env python3
"""Task 3: check of the symmetry-reduction certificate cert8.txt (levels 1, 2).

(a) listed cosets of each node == Cand(W_d, R) = {cosets in C(W_d) meeting R_R}  (from definitions)
(b) each line: U invertible, W_d^U = W_d, U R = R, c_o^U == x (mod W_d)
(c) the 486 spaces W_2 + <c> == the 486 roots of list3_basis.txt (as subspaces)
(d) R = smallest a with dimP=1 if any, else smallest a with dimP=0  (and dimP(R) <= 1)
Action: B^U(x,y) = B(Ux,Uy); (U v)_i = parity(u_i & v); u_i = i-th byte (first hex pair = u_0).
"""
import os
DATA = os.environ.get('N8_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'))
import sys, time
import numpy as np
from indep_lib import (PAIRS, PHI_T, span, echelon, is_good, expand, dimP_direct, bilinear,
                       rank_gf2)

CERT = os.path.join(DATA, 'cert8.txt')
LIST3 = os.path.join(DATA, 'list3_basis.txt')
J = 4202497   # e0^e1 + e2^e3 + e4^e5


def parity(x):
    return (np.bitwise_count(x) & 1).astype(np.int64)


def phi_gen(m, a):
    """phi_a(m) = M_m a, elementwise for arrays m, a."""
    m = np.asarray(m, dtype=np.int64)
    a = np.asarray(a, dtype=np.int64)
    return (PHI_T[0][m & 255, a] ^ PHI_T[1][(m >> 8) & 255, a] ^
            PHI_T[2][(m >> 16) & 255, a] ^ PHI_T[3][(m >> 24) & 255, a]).astype(np.int64)


def phi_pt(m, a):
    return phi_gen(m, np.full(len(m), a, dtype=np.int64))


def cols_of(rows):
    """rows (N,8) of U -> cols (N,8): col_j = U e_j, bit k of col_j = bit j of u_k."""
    rows = rows.astype(np.int64)
    cols = np.zeros_like(rows)
    for j in range(8):
        for k in range(8):
            cols[:, j] |= ((rows[:, k] >> j) & 1) << k
    return cols


def act(m, rows):
    """mask of B^U for each (m[n], U_n): bit (i,j) = B(U e_i, U e_j) = parity(col_i & M col_j)."""
    cols = cols_of(rows)
    ph = [phi_gen(m, cols[:, j]) for j in range(8)]
    out = np.zeros(len(m), dtype=np.int64)
    for p, (i, j) in enumerate(PAIRS):
        out |= parity(cols[:, i] & ph[j]) << p
    return out


def apply_vec(rows, v):
    """U v for scalar v, each U: (N,) ints."""
    rows = rows.astype(np.int64)
    out = np.zeros(len(rows), dtype=np.int64)
    for i in range(8):
        out |= parity(rows[:, i] & v) << i
    return out


def gf2_rank_rows(rows):
    """rank of each 8x8 matrix (rows (N,8)) over F_2, vectorized elimination."""
    R = rows.astype(np.int64).copy()
    N = len(R)
    used = np.zeros((N, 8), dtype=bool)
    rank = np.zeros(N, dtype=np.int64)
    ar = np.arange(N)
    for j in range(8):
        has = (((R >> j) & 1) == 1) & ~used
        anyp = has.any(1)
        piv = np.argmax(has, 1)
        prow = R[ar, piv]
        used[ar[anyp], piv[anyp]] = True
        rank += anyp
        bitj = ((R >> j) & 1) == 1
        elim = bitj & anyp[:, None]
        elim[ar, piv] = False
        R ^= np.where(elim, prow[:, None], 0)
    return rank


def canon_arr(x, vecs, pivs):
    x = np.asarray(x, dtype=np.int64).copy()
    for v, p in zip(vecs, pivs):
        x ^= ((x >> p) & 1) * v
    return x


def cand_set(basis, R):
    """Cand(W,R) from the definitions; returns (sorted canonical reps, |C(W)|, dimP, vecs, pivs)."""
    d = len(basis)
    vecs, pivs = echelon(basis)
    Wels = span(basis)
    dimP = dimP_direct(basis)
    full = [a for a in range(1, 256) if dimP[a] == 2]
    pa = sorted(pivs)
    nrep = 1 << (28 - d)
    out = []
    nC = 0
    CH = 1 << 22
    for s in range(0, nrep, CH):
        z = np.arange(s, min(nrep, s + CH), dtype=np.int64)
        y = expand(z, pa)
        ok = np.ones(len(y), dtype=bool)
        for w in Wels:
            ok &= is_good(y ^ w)
        y = y[ok]
        ok = np.ones(len(y), dtype=bool)
        for a in full:
            for w in Wels:
                ok &= phi_pt(y ^ w, a) != 0
        y = y[ok]
        nC += len(y)
        meet = np.zeros(len(y), dtype=bool)
        for w in Wels:
            meet |= phi_pt(y ^ w, R) == 0
        out.append(y[meet])
    return np.sort(np.concatenate(out)), nC, dimP, vecs, pivs, full


def selftest(rng):
    # act() against the definition B^U(e_i,e_j) = B(U e_i, U e_j) with U v computed from rows
    for _ in range(200):
        while True:
            rows = rng.integers(0, 256, 8)
            if rank_gf2([int(r) for r in rows]) == 8:
                break
        m = int(rng.integers(0, 1 << 28))
        Uv = lambda v: sum((bin(int(rows[i]) & v).count('1') & 1) << i for i in range(8))
        want = 0
        for p, (i, j) in enumerate(PAIRS):
            want |= bilinear(m, Uv(1 << i), Uv(1 << j)) << p
        got = int(act(np.array([m]), rows.reshape(1, 8))[0])
        assert got == want, (m, rows)
        assert gf2_rank_rows(rows.reshape(1, 8))[0] == 8
        # B^U(x,y) = B(Ux,Uy) for random x,y
        x, y = int(rng.integers(1, 256)), int(rng.integers(1, 256))
        assert bilinear(got, x, y) == bilinear(m, Uv(x), Uv(y))
    # rank function against rank_gf2 on random matrices
    rows = rng.integers(0, 256, (2000, 8))
    rk = gf2_rank_rows(rows)
    for n in range(2000):
        assert rk[n] == rank_gf2([int(r) for r in rows[n]])
    print("selftest ok (act vs definition, rank)", flush=True)


def parse_cert():
    nodes = []
    cur = None
    if os.path.exists(CERT):
        f = open(CERT)
    else:                      # the repository ships the certificate compressed
        import lzma
        f = lzma.open(CERT + '.xz', 'rt')
    with f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == 'NODE':
                d = int(t[1])
                basis = [int(v) for v in t[3:3 + d]]
                assert t[3 + d] == 'r' and t[5 + d] == 'reps'
                R = int(t[4 + d])
                K = int(t[6 + d])
                reps = [int(v) for v in t[7 + d:7 + d + K]]
                assert len(reps) == K and len(t) == 7 + d + K
                cur = dict(d=d, hash=t[2], basis=basis, R=R, reps=reps, x=[], o=[], U=[])
                nodes.append(cur)
            else:
                assert len(t) == 3, line
                cur['x'].append(int(t[0]))
                cur['o'].append(int(t[1]))
                cur['U'].append(bytes.fromhex(t[2]))
    for nd in nodes:
        nd['x'] = np.array(nd['x'], dtype=np.int64)
        nd['o'] = np.array(nd['o'], dtype=np.int64)
        nd['U'] = np.frombuffer(b''.join(nd['U']), dtype=np.uint8).reshape(-1, 8).astype(np.int64)
    return nodes


def key_space(basis):
    return tuple(sorted(span(basis)))


def main():
    t0 = time.time()
    rng = np.random.default_rng(7)
    selftest(rng)
    nodes = parse_cert()
    print(f"parsed {len(nodes)} NODE entries, lines per node: {[len(n['x']) for n in nodes]} "
          f"({time.time()-t0:.1f}s)", flush=True)
    allok = True
    for ni, nd in enumerate(nodes):
        d, B, R, reps = nd['d'], nd['basis'], nd['R'], nd['reps']
        tn = time.time()
        print(f"\n== node {ni}: d={d} basis={B} R={R} reps={len(reps)} lines={len(nd['x'])}", flush=True)
        assert rank_gf2(B) == d
        allgood = bool(is_good(np.array([e for e in span(B) if e])).all())
        cand, nC, dimP, vecs, pivs, full = cand_set(B, R)
        # (d)
        ones = [a for a in range(1, 256) if dimP[a] == 1]
        zeros = [a for a in range(1, 256) if dimP[a] == 0]
        Rrule = ones[0] if ones else zeros[0]
        print(f"  all nonzero elements rank 6: {allgood}; dimP histogram {np.bincount(dimP[1:], minlength=3).tolist()}; "
              f"full points {full}")
        print(f"  (d) dimP(R)={dimP[R]}  rule gives R={Rrule}  -> {'OK' if Rrule == R else 'MISMATCH'}")
        # (a)
        xc = canon_arr(nd['x'], vecs, pivs)
        xs = np.sort(xc)
        ndup = int((np.diff(xs) == 0).sum())
        inW = int((xc == 0).sum())
        missing = np.setdiff1d(cand, xs)
        extra = np.setdiff1d(xs, cand)
        okA = ndup == 0 and inW == 0 and len(missing) == 0 and len(extra) == 0 and len(xs) == len(cand)
        print(f"  (a) |C(W)|={nC}  |Cand(W,R)|={len(cand)}  listed={len(xs)}  duplicates={ndup}  "
              f"in-W={inW}  missing={len(missing)}  extra={len(extra)}  -> {'OK' if okA else 'FAIL'}")
        # reps: distinct cosets, each in Cand, own line maps to itself
        rc = canon_arr(np.array(reps), vecs, pivs)
        repsdistinct = len(np.unique(rc)) == len(rc)
        repsincand = bool(np.isin(rc, cand).all())
        print(f"  reps distinct mod W: {repsdistinct}; reps in Cand: {repsincand}; "
              f"o range ok: {bool(((nd['o'] >= 0) & (nd['o'] < len(reps))).all())}; "
              f"orbit sizes: min {np.bincount(nd['o']).min()} max {np.bincount(nd['o']).max()}")
        # (b)
        U = nd['U']
        rk = gf2_rank_rows(U)
        okInv = bool((rk == 8).all())
        Wset = np.array(sorted(span(B)), dtype=np.int64)
        okW = np.ones(len(U), dtype=bool)
        for b in B:
            img = act(np.full(len(U), b, dtype=np.int64), U)
            okW &= np.isin(img, Wset)
        UR = apply_vec(U, R)
        okR = UR == R
        co = np.array(reps, dtype=np.int64)[nd['o']]
        img = act(co, U)
        okX = canon_arr(img ^ nd['x'], vecs, pivs) == 0
        okB = okInv and okW.all() and okR.all() and okX.all()
        print(f"  (b) U invertible: {int((rk==8).sum())}/{len(U)}; W^U=W: {int(okW.sum())}/{len(U)}; "
              f"UR=R: {int(okR.sum())}/{len(U)}; c_o^U = x mod W: {int(okX.sum())}/{len(U)}  "
              f"-> {'OK' if okB else 'FAIL'}   ({time.time()-tn:.1f}s)", flush=True)
        nd['okA'], nd['okB'], nd['okD'] = okA, okB, (Rrule == R and dimP[R] <= 1)
        allok &= okA and okB
    # structure: node 0 is d=1 <J>, nodes 1..3 are W1 + <c> for the 3 reps
    n1 = nodes[0]
    assert n1['d'] == 1
    print(f"\nW1 = <J>: {key_space(n1['basis']) == key_space([J])}")
    want2 = sorted(key_space(n1['basis'] + [c]) for c in n1['reps'])
    have2 = sorted(key_space(n['basis']) for n in nodes[1:])
    print(f"d=2 nodes == W1 + <c> for the d=1 reps (as subspaces): {want2 == have2}")
    # (c)
    sp = []
    for n in nodes[1:]:
        for c in n['reps']:
            sp.append(key_space(n['basis'] + [c]))
    roots = []
    with open(LIST3) as f:
        for line in f:
            t = line.split()
            if t and t[0] == 'L3':
                assert t[4] == 'basis'
                bb = [int(v) for v in t[5:8]]
                assert rank_gf2(bb) == 3
                roots.append((int(t[1]), key_space(bb)))
    rootkeys = [k for _, k in roots]
    okC = (len(sp) == 486 and len(set(sp)) == 486 and len(rootkeys) == 486 and
           len(set(rootkeys)) == 486 and set(sp) == set(rootkeys))
    print(f"(c) spaces W2+<c>: {len(sp)} (distinct {len(set(sp))}); roots in list3: {len(rootkeys)} "
          f"(distinct {len(set(rootkeys))}); idx 0..485: {sorted(i for i,_ in roots) == list(range(486))}; "
          f"equal as sets of subspaces: {set(sp) == set(rootkeys)}  -> {'OK' if okC else 'FAIL'}")
    # all roots consist of rank-6 forms
    print(f"all 486 roots: every nonzero element rank 6: "
          f"{all(bool(is_good(np.array([e for e in k if e])).all()) for k in rootkeys)}")
    print(f"\nSUMMARY (a) {[n['okA'] for n in nodes]}  (b) {[n['okB'] for n in nodes]}  (c) {okC}  "
          f"(d) {[n['okD'] for n in nodes]}   total {time.time()-t0:.1f}s")
    ok = okC and all(bool(n['okA']) and bool(n['okB']) and bool(n['okD']) for n in nodes)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
