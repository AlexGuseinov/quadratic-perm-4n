import numpy as np, time, sys
from indep_lib import *
roots = {0:[4194965, 8852, 1085440], 5:[4194965, 8852, 1245205], 17:[4194965, 8852, 8589332]}
for rid, B in roots.items():
    t=time.time()
    els = span(B)
    print(rid, "rank of basis", rank_gf2(B), "all nonzero good:", bool(is_good(np.array([e for e in els if e])).all()))
    vecs, pivs = echelon(B)
    Y, I, dimP = brute_C(B, pivs)
    print(" pivs", pivs, "dimP hist", np.bincount(dimP[1:]), "|C|", len(Y), f"{time.time()-t:.1f}s")
    cnt = unpack256(I).sum(0)
    cand = np.nonzero(dimP[1:]<=1)[0]+1
    key = 4*cnt[cand] + (dimP[cand]==0)
    r = cand[np.argmin(key)]
    print(" r", r, "cnt[r]", cnt[r], "min cnt", cnt[cand].min(), "max", cnt[cand].max())
