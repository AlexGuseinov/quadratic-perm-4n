"""Same D(n) model as du4perm_sat.py, solved with CryptoMiniSat using native XOR
constraints (Gauss-Jordan elimination) instead of Tseitin-encoded XOR gates.

Usage: python3 du4perm_cms.py n [--sb3] [--noperm] [--threads T] [--time S]
"""
import sys
import time

import pycryptosat

import os, sys  # repository layout: shared modules are in ../../core
sys.path.insert(1, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "core"))
import du4perm_sat as base
from du4perm_sat import T, F_


class XorEnc(base.Enc):
    """Records o = a xor b as a native XOR constraint."""

    def __init__(self):
        super().__init__()
        self.xors = []  # (vars, rhs)

    def xor2(self, a, b):
        if a is F_:
            return b
        if b is F_:
            return a
        if a is T:
            return self.neg(b)
        if b is T:
            return self.neg(a)
        o = self.new()
        # o + a + b = 0  (with negations folded into rhs)
        rhs = False
        vs = []
        for lit in (o, a, b):
            if lit < 0:
                rhs = not rhs
            vs.append(abs(lit))
        self.xors.append((vs, rhs))
        return o


def main():
    n = int(sys.argv[1])
    sb3 = "--sb3" in sys.argv
    perm = "--noperm" not in sys.argv
    threads = int(sys.argv[sys.argv.index("--threads") + 1]) if "--threads" in sys.argv else 1
    tlim = float(sys.argv[sys.argv.index("--time") + 1]) if "--time" in sys.argv else None
    base.Enc = XorEnc  # build() instantiates Enc
    t0 = time.time()
    enc, c, ell = base.build(n, sb=True, sb3=sb3, perm=perm)
    print(f"n={n} sb3={sb3} perm={perm}: vars={enc.pool.top} clauses={len(enc.clauses)} "
          f"xors={len(enc.xors)} (build {time.time() - t0:.1f}s)", flush=True)
    kw = {"threads": threads}
    if tlim:
        kw["time_limit"] = tlim
    s = pycryptosat.Solver(**kw)
    for cl in enc.clauses:
        s.add_clause(cl)
    for vs, rhs in enc.xors:
        s.add_xor_clause(vs, rhs)
    t1 = time.time()
    sat, sol = s.solve()
    dt = time.time() - t1
    if sat is True:
        model = [i if sol[i] else -i for i in range(1, len(sol))]
        C, L = base.extract(model, c, ell, n)
        table = base.lut(C, L, n)
        p, du = base.check(table, n)
        print(f"SAT in {dt:.1f}s: permutation={p} DU={du}")
        print("LUT:", table)
    elif sat is False:
        print(f"UNSAT in {dt:.1f}s")
    else:
        print(f"UNKNOWN after {dt:.1f}s")


if __name__ == "__main__":
    main()
