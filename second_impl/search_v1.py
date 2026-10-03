#!/usr/bin/env python3
"""Task 2: node count of the search tree by the second implementation (numpy implementation).

Usage: python3 search.py ROOTID [--verify] [--timelimit SECONDS]

Data kept per node W (dim d):
  basis  : list of masks spanning W (root basis + one coset rep per level)
  pivs   : pivot positions; the canonical representative of a coset y+W is the unique
           element that is 0 on all pivot positions
  S      : canonical representatives of the cosets in C(W)            (uint32, N)
  I      : incidence bitsets, I[n] = { a : coset S[n]+W meets R_a }    ((N,4) uint64)
           (R_a = forms whose radical contains a; "meets" = some element of the coset
           has a in its radical)
  dimP   : dimP_W(a) for all points a (index 0 unused)

Root: everything computed from the definitions (brute force over all 2^25 cosets).
Child W' = W + <c>, c+W in C(W), new pivot p = highest set bit of c:
  dimP_{W'}(a) = dimP_W(a) + [c+W meets R_a]          (P_a(W') = P_a(W) u ((c+W) n R_a))
  C(W') = { y+W' : y+W in C(W), y+c+W in C(W) } minus cosets meeting R_a for the points a
          that are full for W' but not for W; incidence of y+W' = I(y) | I(y+c).
  canonical rep of y+W' = the one of y, y^c with bit p = 0.
"""
import sys, time, argparse
import numpy as np
from indep_lib import (span, rank_gf2, echelon, is_good, brute_C, dimP_direct, unpack256,
                       points_to_mask, compress)

ROOTS = {0: [4194965, 8852, 1085440],
         5: [4194965, 8852, 1245205],
         17: [4194965, 8852, 8589332]}

THR = {d: (1 << (8 - d)) - 1 for d in range(3, 8)}
POS = {}
LIMIT = 1 << 22   # max elements of a (children x N) work matrix


def pos_table(d):
    if d not in POS:
        POS[d] = np.full(1 << (28 - d), -1, dtype=np.int32)
    return POS[d]


class TimeUp(Exception):
    pass


class Stats:
    def __init__(self):
        z = lambda: [0] * 9
        self.nodes, self.pdim, self.psize, self.pcov, self.expanded = z(), z(), z(), z(), z()
        self.created, self.processed, self.rule3 = z(), z(), z()
        self.solutions = []
        self.root_children_done = 0
        self.root_children_total = 0


ST = Stats()
ARGS = None
T0 = None
VERIFY_CREATED = {4: {1, 2000, 12000}, 5: {1, 30000, 400000, 1000000}, 6: {1, 300, 5000}, 7: {1, 40}}
VERIFY_PROCESSED = {5: {1, 50, 2000}, 6: {1, 20}, 7: {1}}
VERIFY_RULE3 = {5: {1, 2, 30}, 6: {1, 10}, 7: {1}}
VLOG = []


def log(msg):
    print(msg, flush=True)


def colsum(I):
    tot = np.zeros(256, np.int64)
    for s in range(0, len(I), 1 << 15):
        tot += unpack256(I[s:s + (1 << 15)]).sum(0, dtype=np.int64)
    return tot


def verify_node(tag, basis, pivs, S, I, dimP):
    """Compare incrementally obtained (S, I, dimP) with a from-scratch computation."""
    t = time.time()
    Sb, Ib, dPb = brute_C(basis, pivs)
    o = np.argsort(S)
    okS = len(Sb) == len(S) and np.array_equal(Sb.astype(np.int64), S[o].astype(np.int64))
    okI = okS and np.array_equal(Ib, I[o])
    okD = np.array_equal(dPb[1:], dimP[1:])
    msg = (f"VERIFY {tag} d={len(basis)} |C| inc={len(S)} brute={len(Sb)} "
           f"S:{okS} I:{okI} dimP:{okD} ({time.time()-t:.1f}s)")
    VLOG.append((okS and okI and okD, msg))
    log(msg)
    if not (okS and okI and okD):
        log("   basis=" + str(basis) + " pivs=" + str(pivs))


def process(d, basis, pivs, S, I, dimP):
    """Process a node of dimension d (already counted in ST.nodes[d])."""
    if d == 8:
        if (dimP[1:] == 2).all():
            ST.solutions.append(list(basis))
            dP = dimP_direct(basis)
            log(f"SOLUTION {basis} direct-dimP-all-2={bool((dP[1:]==2).all())} "
                f"allgood={bool(is_good(np.array([e for e in span(basis) if e])).all())}")
        else:
            ST.pdim[8] += 1
        return
    ST.processed[d] += 1
    if ARGS.verify and ST.processed[d] in VERIFY_PROCESSED.get(d, ()):
        verify_node(f"processed#{ST.processed[d]}", basis, pivs, S, I, dimP)
    N = len(S)
    if N < THR[d]:
        ST.psize[d] += 1
        return
    ST.rule3[d] += 1
    if ARGS.verify and ST.rule3[d] in VERIFY_RULE3.get(d, ()):
        verify_node(f"reached-rule3#{ST.rule3[d]}", basis, pivs, S, I, dimP)
    cnt = colsum(I)
    cand = np.nonzero(dimP[1:] <= 1)[0] + 1          # increasing a
    cc = cnt[cand]
    dz = (dimP[cand] == 0)
    if (cc == 0).any() or (dz & (cc < 3)).any():
        ST.pcov[d] += 1
        return
    key = 4 * cc + dz.astype(np.int64)
    r = int(cand[np.argmin(key)])                     # first minimum = smallest a
    ST.expanded[d] += 1
    rows = np.nonzero((I[:, r >> 6] >> np.uint64(r & 63)) & np.uint64(1))[0]
    K = len(rows)
    assert K == cnt[r]
    ST.nodes[d + 1] += K
    if d == 3:
        ST.root_children_total = K
    cdim = dimP[None, :] + unpack256(I[rows]).astype(np.int64)
    cdim[:, 0] = 0
    isdim = (cdim[:, 1:] >= 3).any(1)
    ST.pdim[d + 1] += int(isdim.sum())
    C = S[rows]
    if d + 1 == 8:
        for k in range(K):
            if not isdim[k]:
                process(8, basis + [int(C[k])], None, None, None, cdim[k])
        return
    live = np.nonzero(~isdim)[0]
    T = pos_table(d)
    pd = sorted(pivs, reverse=True)
    cS = compress(S, pd)
    T[cS] = np.arange(N, dtype=np.int32)
    thr = THR[d + 1]
    Kc = max(1, LIMIT // max(N, 1))
    try:
        for s in range(0, len(live), Kc):
            ks = live[s:s + Kc]
            Cc = C[ks]
            P = (np.frexp(Cc.astype(np.float64))[1] - 1).astype(np.uint32)
            Y = S[None, :] ^ Cc[:, None]
            sel = ((S[None, :] >> P[:, None]) & np.uint32(1)) == 0
            idx = T[compress(Y, pd)]
            valid = sel & (idx >= 0)
            del Y, sel
            npairs = valid.sum(1)
            for t, k in enumerate(ks):
                ST.created[d + 1] += 1
                c = int(Cc[t]); p = int(P[t])
                verify_this = ARGS.verify and ST.created[d + 1] in VERIFY_CREATED.get(d + 1, ())
                if npairs[t] < thr and not verify_this:
                    ST.psize[d + 1] += 1
                else:
                    ii = np.nonzero(valid[t])[0]
                    jj = idx[t, ii]
                    nS = S[ii]
                    nI = I[ii] | I[jj]
                    newfull = (cdim[k] == 2) & (dimP < 2)
                    newfull[0] = False
                    if newfull.any():
                        F = points_to_mask(newfull)
                        keep = ~((nI & F) != 0).any(1)
                        nS = nS[keep]
                        nI = nI[keep]
                    if verify_this:
                        verify_node(f"created#{ST.created[d+1]}", basis + [c], pivs + [p], nS, nI, cdim[k])
                    if npairs[t] < thr:
                        ST.psize[d + 1] += 1
                    else:
                        process(d + 1, basis + [c], pivs + [p], nS, nI, cdim[k])
                if d == 3:
                    ST.root_children_done += 1
                    if ST.root_children_done % 500 == 0:
                        report(partial=True)
                    if ARGS.timelimit and time.time() - T0 > ARGS.timelimit:
                        raise TimeUp()
    finally:
        T[cS] = -1


def report(partial):
    el = time.time() - T0
    hdr = "PARTIAL" if partial else "FINAL"
    lines = [f"[{hdr}] root-children done {ST.root_children_done}/{ST.root_children_total}  elapsed {el:.1f}s"]
    lines.append("level   nodes   pruned-dim pruned-size pruned-cov  expanded  (processed)")
    for d in range(3, 9):
        lines.append(f"  {d}  {ST.nodes[d]:9d} {ST.pdim[d]:10d} {ST.psize[d]:11d} {ST.pcov[d]:10d} "
                     f"{ST.expanded[d]:9d}  ({ST.processed[d]})")
    lines.append(f"  solutions: {len(ST.solutions)}")
    log("\n".join(lines))


def main():
    global ARGS, T0
    ap = argparse.ArgumentParser()
    ap.add_argument('root', type=int)
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--timelimit', type=float, default=0)
    ARGS = ap.parse_args()
    T0 = time.time()
    B = ROOTS[ARGS.root]
    els = span(B)
    assert rank_gf2(B) == 3, "root basis not independent"
    allgood = bool(is_good(np.array([e for e in els if e])).all())
    log(f"root {ARGS.root}: basis {B}  dim={rank_gf2(B)}  all 7 nonzero elements rank 6: {allgood}")
    assert allgood
    vecs, pivs = echelon(B)
    S, I, dimP = brute_C(B, pivs)          # from the definition
    S = S.astype(np.uint32)
    log(f"root: pivots {pivs}  dimP histogram {np.bincount(dimP[1:], minlength=4).tolist()}  "
        f"|C(W)| = {len(S)}  ({time.time()-T0:.1f}s)")
    ST.nodes[3] = 1
    if (dimP[1:] >= 3).any():
        ST.pdim[3] += 1
    else:
        try:
            process(3, list(B), list(pivs), S, I, dimP)
            report(partial=False)
        except TimeUp:
            log("TIME LIMIT REACHED -- counts below are PARTIAL")
            report(partial=True)
    for s in ST.solutions:
        log(f"solution basis: {s}")
    if VLOG:
        log(f"verifications: {sum(ok for ok,_ in VLOG)}/{len(VLOG)} passed")
    log(f"total time {time.time()-T0:.1f}s")


if __name__ == '__main__':
    main()
