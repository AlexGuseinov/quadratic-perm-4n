#!/usr/bin/env python3
"""Task 2/4: node count of the search tree by the second implementation (numpy implementation).

Single root : python3 search.py ROOTID [--verify] [--timelimit SECONDS]
              (ROOTID = idx in list3_basis.txt)
Batch       : python3 search.py --batch ROOTS_FILE --out RESULTS_FILE [--start I] [--step M]
              processes the roots with idx % M == I, appends one line per finished root,
              skips roots whose idx already has a 'ROOT idx=...' line in RESULTS_FILE.

Data kept per node W (dim d):
  basis  : list of masks spanning W (root basis + one coset rep per level)
  pivs   : pivot positions; the canonical representative of a coset y+W is the unique
           element that is 0 on all pivot positions
  S      : SORTED canonical representatives of the cosets in C(W)     (uint32, N)
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
Membership "y^c in C(W)" is tested with a per-level bitmap over all 2^28 masks.
(search_v1.py is the earlier version -- pos tables instead of bitmaps -- that produced
 root0.log, root5.log, root17.log; the counts are identical.)
"""
import os
DATA = os.environ.get('N8_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data'))
import os, sys, time, argparse
import numpy as np
from indep_lib import (span, rank_gf2, echelon, is_good, brute_C, dimP_direct, unpack256,
                       points_to_mask)

ROOTS_DEFAULT = os.path.join(DATA, 'list3_basis.txt')
THR = {d: (1 << (8 - d)) - 1 for d in range(3, 8)}
LIMIT = 1 << 22   # max elements of a (children x N) work matrix
BM = {}           # per-level byte bitmaps over all 2^28 masks (32 MB each, levels 3..7)


def bitmap(d):
    if d not in BM:
        BM[d] = np.zeros(1 << 25, dtype=np.uint8)
    return BM[d]


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
QUIET = False
VERIFY_CREATED = {4: {1, 2000, 12000}, 5: {1, 30000, 400000, 1000000}, 6: {1, 300, 5000}, 7: {1, 40}}
VERIFY_PROCESSED = {5: {1, 50, 2000}, 6: {1, 20}, 7: {1}}
VERIFY_RULE3 = {5: {1, 2, 30}, 6: {1, 10}, 7: {1}}
VLOG = []


def log(msg, force=False):
    if force or not QUIET:
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
    log(msg, force=True)
    if not (okS and okI and okD):
        log("   basis=" + str(basis) + " pivs=" + str(pivs), force=True)


def process(d, basis, pivs, S, I, dimP):
    """Process a node of dimension d (already counted in ST.nodes[d])."""
    if d == 8:
        if (dimP[1:] == 2).all():
            ST.solutions.append(list(basis))
            dP = dimP_direct(basis)
            log(f"SOLUTION {basis} direct-dimP-all-2={bool((dP[1:]==2).all())} "
                f"allgood={bool(is_good(np.array([e for e in span(basis) if e])).all())}", force=True)
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
    P_all = (np.frexp(C.astype(np.float64))[1] - 1).astype(np.int64)   # highest bit of c
    B8 = bitmap(d)
    np.bitwise_or.at(B8, S >> 3, (np.uint8(1) << (S & 7).astype(np.uint8)))
    thr = THR[d + 1]
    try:
        for p in np.unique(P_all[live]):
            p = int(p)
            rows0 = np.nonzero(((S >> p) & 1) == 0)[0]     # y with bit p = 0
            S0 = S[rows0]
            N0 = max(len(S0), 1)
            ksp = live[P_all[live] == p]
            Kc = max(1, LIMIT // N0)
            for s in range(0, len(ksp), Kc):
                ks = ksp[s:s + Kc]
                Cc = C[ks]
                Y = S0[None, :] ^ Cc[:, None]
                hit = ((B8[Y >> 3] >> (Y & 7).astype(np.uint8)) & 1).astype(bool)
                npairs = hit.sum(1)
                c0 = ST.created[d + 1]
                ST.created[d + 1] += len(ks)
                vt = set()
                if ARGS.verify:
                    vt = {v - c0 - 1 for v in VERIFY_CREATED.get(d + 1, ()) if c0 < v <= c0 + len(ks)}
                work = np.nonzero(npairs >= thr)[0].tolist()
                ST.psize[d + 1] += int((npairs < thr).sum())
                for t in sorted(set(work) | vt):
                    k = ks[t]
                    c = int(Cc[t])
                    ii0 = np.nonzero(hit[t])[0]
                    ii = rows0[ii0]
                    jj = np.searchsorted(S, Y[t, ii0])
                    assert np.array_equal(S[jj], Y[t, ii0])
                    nS = S[ii]
                    nI = I[ii] | I[jj]
                    newfull = (cdim[k] == 2) & (dimP < 2)
                    newfull[0] = False
                    if newfull.any():
                        F = points_to_mask(newfull)
                        keep = ~((nI & F) != 0).any(1)
                        nS = nS[keep]
                        nI = nI[keep]
                    if t in vt:
                        verify_node(f"created#{c0+t+1}", basis + [c], pivs + [p], nS, nI, cdim[k])
                    if npairs[t] >= thr:
                        process(d + 1, basis + [c], pivs + [p], nS, nI, cdim[k])
                if d == 3:
                    before = ST.root_children_done
                    ST.root_children_done += len(ks)
                    if ST.root_children_done // 2000 != before // 2000:
                        report(partial=True)
                    if ARGS.timelimit and time.time() - T0 > ARGS.timelimit:
                        raise TimeUp()
    finally:
        B8[S >> 3] = 0


def report(partial):
    el = time.time() - T0
    hdr = "PARTIAL" if partial else "FINAL"
    lines = [f"[{hdr}] root-children done {ST.root_children_done}/{ST.root_children_total}  elapsed {el:.1f}s"]
    lines.append("level   nodes   pruned-dim pruned-size pruned-cov  expanded  (processed)")
    for d in range(3, 9):
        lines.append(f"  {d}  {ST.nodes[d]:9d} {ST.pdim[d]:10d} {ST.psize[d]:11d} {ST.pcov[d]:10d} "
                     f"{ST.expanded[d]:9d}  ({ST.processed[d]})")
    lines.append(f"  solutions: {len(ST.solutions)}")
    log("\n".join(lines), force=not partial)


def run_root(B, label):
    """Run the search from root basis B. Returns True if complete, False if time limit hit."""
    global ST, T0
    ST = Stats()
    T0 = time.time()
    els = span(B)
    assert rank_gf2(B) == 3, "root basis not independent"
    allgood = bool(is_good(np.array([e for e in els if e])).all())
    log(f"root {label}: basis {B}  dim={rank_gf2(B)}  all 7 nonzero elements rank 6: {allgood}", force=True)
    assert allgood
    vecs, pivs = echelon(B)
    S, I, dimP = brute_C(B, pivs)          # from the definition (sorted)
    S = S.astype(np.uint32)
    log(f"root: pivots {pivs}  dimP histogram {np.bincount(dimP[1:], minlength=4).tolist()}  "
        f"|C(W)| = {len(S)}  ({time.time()-T0:.1f}s)", force=True)
    ST.nodes[3] = 1
    if (dimP[1:] >= 3).any():
        ST.pdim[3] += 1
        return True
    try:
        process(3, list(B), list(pivs), S, I, dimP)
        return True
    except TimeUp:
        log("TIME LIMIT REACHED -- counts are PARTIAL", force=True)
        return False


def read_roots(fn):
    roots = []
    with open(fn) as f:
        for line in f:
            t = line.split()
            if t and t[0] == 'L3':
                assert t[4] == 'basis', line
                roots.append((int(t[1]), [int(v) for v in t[5:8]]))
    return roots


def result_line(tag, idx, B, secs):
    L = range(3, 9)
    f = lambda arr: ",".join(str(arr[d]) for d in L)
    s = (f"{tag} idx={idx} basis={','.join(map(str, B))} nodes={f(ST.nodes)} pdim={f(ST.pdim)} "
         f"psize={f(ST.psize)} pcov={f(ST.pcov)} expanded={f(ST.expanded)} "
         f"sols={len(ST.solutions)} secs={secs:.1f}")
    if ST.solutions:
        s += " solbases=" + ";".join(":".join(map(str, sb)) for sb in ST.solutions)
    return s


def done_set(fn):
    done = set()
    if os.path.exists(fn):
        with open(fn) as f:
            for line in f:
                if line.startswith("ROOT idx="):
                    done.add(int(line.split()[1][4:]))
    return done


def batch():
    global QUIET
    QUIET = True
    roots = read_roots(ARGS.batch)
    mine = [(i, B) for i, B in roots if i % ARGS.step == ARGS.start]
    done = done_set(ARGS.out)
    todo = [(i, B) for i, B in mine if i not in done]
    log(f"batch: {len(roots)} roots in file, {len(mine)} with idx % {ARGS.step} == {ARGS.start}, "
        f"{len(mine)-len(todo)} already in {ARGS.out}, {len(todo)} to do", force=True)
    tb = time.time()
    for n, (i, B) in enumerate(todo):
        if ARGS.maxroots and n >= ARGS.maxroots:
            break
        t = time.time()
        complete = run_root(B, f"idx={i}")
        secs = time.time() - t
        line = result_line("ROOT" if complete else "PARTIAL", i, B, secs)
        with open(ARGS.out, 'a') as f:        # one short write per line (O_APPEND)
            f.write(line + "\n")
            f.flush()
            os.fsync(f.fileno())
        log(line, force=True)
        log(f"batch progress: {n+1}/{len(todo)} done, elapsed {time.time()-tb:.0f}s", force=True)


def main():
    global ARGS
    ap = argparse.ArgumentParser()
    ap.add_argument('root', type=int, nargs='?')
    ap.add_argument('--roots', default=ROOTS_DEFAULT, help='root list for single-root mode')
    ap.add_argument('--batch', metavar='ROOTS_FILE')
    ap.add_argument('--out', metavar='RESULTS_FILE')
    ap.add_argument('--start', type=int, default=0)
    ap.add_argument('--step', type=int, default=1)
    ap.add_argument('--maxroots', type=int, default=0, help='(testing) stop after this many roots')
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--timelimit', type=float, default=0, help='per-root time limit (s), 0 = none')
    ARGS = ap.parse_args()
    if ARGS.batch:
        assert ARGS.out, "--out required with --batch"
        batch()
        return
    roots = dict(read_roots(ARGS.roots))
    B = roots[ARGS.root]
    complete = run_root(B, ARGS.root)
    report(partial=not complete)
    for s in ST.solutions:
        log(f"solution basis: {s}", force=True)
    if VLOG:
        log(f"verifications: {sum(ok for ok,_ in VLOG)}/{len(VLOG)} passed", force=True)
    log(result_line("ROOT" if complete else "PARTIAL", ARGS.root, B, time.time() - T0), force=True)
    log(f"total time {time.time()-T0:.1f}s", force=True)


if __name__ == '__main__':
    main()
