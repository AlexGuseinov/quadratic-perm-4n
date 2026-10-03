import numpy as np, time
from indep_lib import *
rng = np.random.default_rng(12345)
t=time.time()
# 1. phi-based radical vs definition of B, on random masks
for m in rng.integers(0, 1<<28, 300):
    m=int(m)
    Z = np.nonzero(phi_all(np.array([m]))[0]==0)[0]
    Z = [int(a) for a in Z if a]
    assert Z == radical_points_direct(m), m
print("phi radical == definition radical on 300 random masks")
# 2. Pfaffian table vs radical size (rank = 8 - dim rad) on 2M random masks
M = rng.integers(0, 1<<28, 2_000_000)
g = is_good(M)
nrad = np.zeros(len(M), dtype=np.int64)
for s in range(0, len(M), 1<<16):
    nrad[s:s+(1<<16)] = (phi_all(M[s:s+(1<<16)])[:,1:]==0).sum(1)
r6 = nrad == 3
print("rank6 by Pfaffian vs radical size: mismatches", int((g!=r6).sum()), "of", len(M), " frac rank6", g.mean())
# also exhaustive for masks < 2^20 (all forms supported on low pairs)
M = np.arange(1<<20)
g = is_good(M)
nrad = np.zeros(len(M), dtype=np.int64)
for s in range(0, len(M), 1<<16):
    nrad[s:s+(1<<16)] = (phi_all(M[s:s+(1<<16)])[:,1:]==0).sum(1)
print("exhaustive masks <2^20 mismatches:", int((g!=(nrad==3)).sum()))
print(f"{time.time()-t:.1f}s")
