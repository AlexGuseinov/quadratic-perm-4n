# quadratic-perm-4n

Research code for the question: **does F_2^n with 4 | n admit a quadratic
permutation with differential uniformity 4?** Answer for n = 8: **no** (computer-assisted,
2026-09-28; see below).

Status (2026-09-28): the proofs in `notes/NOTES_lemmas.md` have not yet been checked by the author.

To verify the n = 8 theorem, start with [`ARTIFACT.md`](ARTIFACT.md).

## Results so far

| Result | Status |
|---|---|
| L1 count identity; L2: a DU-4 quadratic permutation has all components of rank n-2 and all derivative kernels of dimension 2 | known: Charpin–Peng, FFA 56 (2019), Thms 5 and 7 (short proof included) |
| L3/L4: derivative kernels form a line spread, totally isotropic for all components, mapped onto a spread | proved (short) |
| L7/L8: the kernel spread cannot be Desarguesian when 4 \| n; general t: Desarguesian t-spread forces n/t odd (t = 1 gives the classical "quadratic APN permutations need n odd") | proved (paper, Theorem `thm:general-t`) |
| L6: F_4-quadratic maps (incl. Gold-type) are never DU-4 permutations when 4 \| n | proved (special case of L7) |
| n = 4: no DU-4 quadratic permutation | MacWilliams proof + SAT with verified DRAT proof |
| n = 8, Desarguesian spread | SAT: UNSAT in 1 s, DRAT proof verified (`earlier_approaches/sat/proofs/d8_desarg.drat`) |
| n = 8: the quadratic part beta of a DU-4 permutation has **no symmetry** beta(Ax,Ay) = B beta(x,y), (A,B) != (I,I) | proved (computer-assisted): 93 pairs (A of prime order up to conjugacy, B up to conjugacy); 79 excluded by Lemma sym (paper, NOTES N1-N10); 14 by kissat with DRAT or exact enumeration, see `logs/sym_summary_n8.txt` |
| n = 8: no two kernel lines span a 4-space that is a union of kernel lines (no closed pair); for every even n >= 6: k(L,M) in {2,3,5} and dim span beta(L,M) = 4, 3, 2 accordingly | proved (paper Prop. prop:local; NOTES K4-K7; `dimension6/check_closed_lemma.py`) |
| **n = 8: no quadratic DU-4 permutation exists** (Alt(8, F_2) has no 8-dim subspace of constant rank 6 in which every point lies in exactly 3 radicals; so every quadratic (8,8)-function with all derivative kernels of dim 2 has a bent component) | **proved (computer-assisted)**: radical-closure search `search_n8/crsub_enum.c` (paper Sec. 5, NOTES "n = 8 by exhaustive search"): 486 nodes of dim 3, 1,466,318 of dim 4, 244,686,904 of dim 5, none of dim 6; about 4 CPU-hours. Rerun from archived sources with identical counts (`search_n8/run_n8.sh`, `logs/crsub_n8_final_*`). Cross-checked by a second, separately written Python/numpy implementation (`second_impl/`: Pfaffian rank table, certificate of the symmetry reduction at levels 1-2, all 486 subtrees searched without symmetry: 0 solutions, `second_impl/batch_summary.txt`) and at n = 6 (finds exactly the two known classes; rule-free enumeration `search_n8/crsub_brute2.c` gives the same 12,288 subspaces containing J) |
| n = 6: Alt(6, F_2) has exactly two GL(6,2)-classes of 6-dim subspaces of constant rank 4; condition (ii) (every point in exactly 3 radicals) follows from constant rank | computed: `search_n8/crsub_brute3.c` (no condition on the P_a) finds the same 12,288 subspaces containing J as `crsub_brute2.c`, all with dim P_a = 2 (`logs/crsub_brute3_n6.log`) |
| n = 8, earlier approaches (SAT on the case split B1/B2 by local structure) | superseded by the search; kissat runs of 1-4 h gave no verdict (`earlier_approaches/sat/n8_cases.py`, logs) |

## Repository layout

```
README.md, ARTIFACT.md       this file; how to verify the n = 8 theorem
expected_counts.json         expected node counts of both programs (machine-readable)
check_all.sh                 quick check of the small computations (5 seconds)
paper/                       the paper (main.tex, main.pdf)
search_n8/                   first program: the n = 8 search (C), n = 6 enumerations
second_impl/                 second program (Python/NumPy), certificate, its logs
core/                        shared Python modules and the quick checks
dimension6/                  experiments on n = 6 solutions
earlier_approaches/          structured classes, SAT, symmetries (superseded for n = 8)
logs/                        outputs of all runs
notes/                       working notes and internal reviews
```

### `search_n8/` — the n = 8 search (Theorem thm:n8)

| File | What it does |
|---|---|
| `crsub_enum.c`, `crsub_common.h` | **The search**: radical-closure search over constant-rank subspaces with (ii); options `--part`, `--nosym`, `--nosym-from`, `--cert`, `--list3`, `--sample` (see the header) |
| `run_n8.sh`, `check_n8_logs.py` | Full n = 8 run with provenance record (sha256 of sources, binary, rank table); completeness check of the logs (compares with `expected_counts.json`) |
| `crsub_sample.c` | Size estimates of the sets of good cosets (independence heuristic vs. exact counts) |
| `crsub_brute.c`, `crsub_brute2.c` | Rule-free enumerations for n = 6 (`crsub_brute2`: each subspace once, via greedy bases) |
| `crsub_brute3.c` | As `crsub_brute2`, but with constant rank as the only condition (no condition on the P_a); n = 6 |
| `crsub_invariants.py` | GL-invariants of found subspaces (kernel-spread closed pairs, radicals, permutation test) |

### `second_impl/` — second implementation

Python/NumPy, written separately: Pfaffian rank table, certificate checker, subtree search, its
logs, and the certificate of the two top levels (`data/cert8.txt.xz`). See `second_impl/README.md`.

### `core/` — shared modules and quick checks

| File | What it does |
|---|---|
| `quadtools.py` | GF(2^n) arithmetic, ranks of alternating forms, DU, sign of a permutation |
| `du4perm_sat.py` | **Exact SAT model for D(n)** with symmetry breaking; options `--sb3`, `--noperm`, `--desarg` |
| `constrank_sat.py` | SAT model: constant-rank subspaces of Alt(n, F_2) (optional s_a >= 2) |
| `validate_basics.py` | Rank counts in Alt(n, F_2); x^5 on GF(2^6); n = 4 constant-rank bound |
| `macwilliams_alt.py` | MacWilliams test for constant-rank subspaces (alternating-forms scheme) |
| `desarguesian_check.py` | Tests whether a kernel spread is Desarguesian; checks the F_4-bilinear Phi |
| `check_x5.py` | Hand-checkable x^5 examples |

### `dimension6/` — experiments for n = 6

| File | What it does |
|---|---|
| `structure_n6.py`, `structure_perm_n6.py`, `observations_n6.py` | Structure experiments on n = 6 solutions |
| `regulus_switch.py` | Regulus switching of the Desarguesian spread; dim Sigma^perp and its ranks |
| `check_closed_lemma.py` | n = 6 check of the ingredients of the no-closed-pair lemma |

### `earlier_approaches/` — the three approaches of Section 6 (now superseded for n = 8)

| File | What it does |
|---|---|
| `structured/constrank8.c` | Exhaustive searches over structured DO families at n = 8 |
| `sat/export_dimacs.py` | Writes the model as DIMACS (for kissat / proof checking) |
| `sat/du4perm_cms.py` | The exact model for CryptoMiniSat with native XOR constraints |
| `sat/n8_cases.py` | Case split B1 / B2 of the general n = 8 problem (options `--local`, `--d4`, `--noperm`) |
| `sat/cube_sample.py`, `sat/cube_runner.py`, `sat/run_d8_half.py`, `sat/vps/` | Cube-and-conquer helpers (difficulty estimates, checkpointed cubes, server scripts) |
| `sat/ti_spread.py`, `sat/ti_spread_multi.py` | Line spreads totally isotropic for one form / a k-dim space of forms |
| `sat/cnf/`, `sat/proofs/` | Small CNFs, `d8_desarg.cnf` and its verified DRAT proof |
| `symmetry/sym_search.py` | Symmetry cases: conjugacy classes, equivariant space, lemma checks and normal forms (`normal_form`) |
| `symmetry/sym_runner.py` | Runs symmetry cases with kissat (+ DRAT check), `--nf` for the normal form |
| `symmetry/sym_twostage.py`, `symmetry/beta_enum.c` | Exact enumeration of the equivariant space (stage 1) + linear part (stage 2) |
| `symmetry/sym_case44.py`, `symmetry/sym_case32.py`, `symmetry/sym_p2.py` | The split cases (p = 5 fixed-point-free; p = 3 case 32 by GL(3,4) classes; involutions) |
| `symmetry/sym_summary.py` | Table of all 93 pairs and how each is excluded (`logs/sym_summary_n8.txt`) |
| `symmetry/validate_sym_n6.py` | n = 6 validation of the whole symmetric pipeline -> `logs/validate_sym_n6.log` |
| `symmetry/conjH_sym.py` | Constant-rank test (Conjecture H) on symmetric subspaces |
| `bin/` | kissat and drat-trim binaries (Linux x86-64) used by the scripts above |
| `build_mac.sh` | Builds drat-trim and `constrank8` on macOS and checks the Desarguesian DRAT proof |

### Other

| File | What it does |
|---|---|
| `paper/main.tex`, `paper/main.pdf` | Paper (19 pages; single file, bibliography inside) |
| `paper/removed_2026-09-29_section6_appendixA.tex` | Text removed when Section 6 was shortened (symmetric-case theorem with proof, structured classes, SAT model, old Appendix A); not compiled |
| `logs/` | All outputs. `logs/NOTES.md`: a reporting bug in one early log, and the folder reorganisation of 2026-10-03 (the logs were written with the earlier flat layout) |
| `notes/NOTES_lemmas.md` | Lemmas L1–L8, SB1–SB4, observations, history of the computation |
| `notes/REVIEW_*.md` | Internal reviews |

The large n = 8 CNFs are not included; regenerate them with
`python3 earlier_approaches/sat/export_dimacs.py 8 d8_sb3.cnf --sb3` (add `--noperm` / `--desarg`).

## Reproducing the non-existence theorem (n = 8)

Details and expected output: `ARTIFACT.md`. All commands are run from the repository root.

```
sh search_n8/run_n8.sh              # builds crsub_enum, 2 parallel parts, about 4.3 CPU-hours;
                                    # writes logs/crsub_n8_final_{part0,part1,provenance}.*
python3 search_n8/check_n8_logs.py  # all 486 level-3 nodes processed, no node of dim 6, PASS
sh second_impl/run_all.sh           # verification WITHOUT crsub_enum: Pfaffian table,
                                    # certificate check, 486 subtrees, PASS
python3 second_impl/summarize_batch.py   # summary of the reported run: 486/486 roots, PASS
```

Dimension 6 (validation of the method), in `search_n8/` after `gcc -O3 -march=native -o X X.c -lm`
for `X` = `crsub_enum`, `crsub_brute2`, `crsub_brute3`:

```
./crsub_enum 6 --out sol6.txt; ./crsub_enum 6 --nosym --out sol6_nosym.txt
python3 crsub_invariants.py 6 sol6.txt --perm           # the two classes
./crsub_brute2 6 | sort -u | wc -l                      # 12288, same set as --nosym (21 min)
./crsub_brute3 6 | cut -d' ' -f2- | sort -u | wc -l     # 12288 again: constant rank alone (5 min)
```

## Quick check of the small computations

`bash check_all.sh` (about 5 seconds; needs `pip install python-sat`). By hand:

```
python3 core/validate_basics.py
python3 core/macwilliams_alt.py
python3 core/du4perm_sat.py 4 --sb3          # UNSAT
python3 core/du4perm_sat.py 6 --sb3          # SAT, prints a permutation with DU 4
python3 earlier_approaches/sat/export_dimacs.py 8 d8_desarg.cnf --sb3 --desarg
kissat d8_desarg.cnf d8_desarg.drat && drat-trim d8_desarg.cnf d8_desarg.drat   # s VERIFIED
```
kissat: https://github.com/arminbiere/kissat, drat-trim: https://github.com/marijnheule/drat-trim

## Reproducing the symmetric-case computation (n = 8; earlier approach)

In `earlier_approaches/symmetry/`:

```
cc -O3 -o beta_enum beta_enum.c
python3 sym_runner.py 8 600 --nf --proof-max 300 --only 13,17,18,27,39,61,64   # kissat + DRAT
python3 sym_twostage.py 8 --nf 54 61 64 65                  # exact enumeration
python3 sym_case44.py; python3 sym_case32.py 1800           # split cases
python3 sym_p2.py 3600 11 16; python3 sym_p2.py 3600 12 --outnorm
python3 sym_summary.py                                      # expect: open: none
python3 validate_sym_n6.py                                  # n = 6 validation
```
Pairs excluded by the lemma need no run (`sym_summary.py` re-checks them). Total CPU time is
about one hour on one core (kissat and drat-trim binaries expected in `earlier_approaches/bin/`).
