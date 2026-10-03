# Artifact guide

How to check the computational statements of the paper
"Quadratic Differentially 4-Uniform Permutations in Dimensions Divisible by Four:
A Desarguesian Obstruction and Non-Existence in Dimension 8" (`paper/main.pdf`).

The mathematical statements the computation relies on are proved in the paper: the reduction to
conditions (i), (ii) (Section 5), the radical-closure lemma (Lemma `lem:closure`) and the
completeness of the search (Appendix A, Theorem `thm:B-complete`). What the programs have to
establish is one fact: **the search tree of Appendix A has no node of dimension 6.**

## 1. What has to be trusted

| Item | Where it comes from | Checked by |
|---|---|---|
| Table of the forms of rank 6 (149,920,960 of 2^28) | computed by each program itself: Gaussian elimination (C), Pfaffians (Python) | closed formula for the number of alternating forms of rank 6; the two tables are bit-identical |
| Certificate of the symmetry reduction at dimensions 1-2 (`second_impl/data/cert8.txt.xz`, 1,468,672 lines) | written by the first program | `second_impl/verify_cert.py` (second program; shares no code with the writer) |
| The 486 nodes of dimension 3 (`second_impl/data/list3_basis.txt`) | written by the first program | `second_impl/verify_cert.py`, part (c): they are exactly the subspaces U + <c_o> of the certificate |
| No node of dimension 6 below the 486 nodes | search | first program (with symmetry reduction), second program (without) |

Minimal trusted base for Theorem `thm:n8`: the proofs in the paper, Python with NumPy, and the
files `second_impl/{indep_lib,task1_pfaffian,verify_cert,search,summarize_batch}.py`
(about 1,050 lines). The first program is then not needed: its only outputs that are used, the
certificate and the list of the 486 nodes, are checked against the definitions.

## 2. Verification without the first program (one command)

```
sh second_impl/run_all.sh [PROCESSES]        # default 2 processes
```

| Step | Script | Time | Expected output |
|---|---|---|---|
| 1 | `task1_pfaffian.py` | 10 s | `rank-6 forms: 149920960`, equal to the closed formula |
| 2 | `verify_cert.py` | 1 min | `SUMMARY (a) [True, ...] (b) [...] (c) True (d) [...]`, then `PASS` |
| 3 | `search.py --batch` | about 1-2 min per node of dimension 3; 7-13 hours of process time in total, 0.6 GB per process | one `ROOT idx=...` line per node, all with `sols=0` |
| 4 | `summarize_batch.py` | 1 s | `roots: 486 of 486`, node counts equal to `expected_counts.json`, `PASS` |

Step 3 can be interrupted and restarted (finished nodes are skipped). Logs of the reported run:
`second_impl/task1.log`, `verify_cert.log`, `batch0.log`, `batch1.log`, `batch_results_*.txt`,
`batch_summary.txt`.

## 3. The first program

```
sh search_n8/run_n8.sh              # builds crsub_enum, two parallel parts, about 4.3 CPU-hours
python3 search_n8/check_n8_logs.py  # 486 nodes of dimension 3 processed, no node of dimension 6,
                                    # node counts equal to expected_counts.json, PASS
```

The sources are in `search_n8/`. `run_n8.sh` records the compiler, the command lines and the SHA-256 hashes of the sources, the
binary and the rank table in `logs/crsub_n8_final_provenance.txt`; the sources of that run are
archived in `logs/crsub_n8_final_src/` (identical to `search_n8/crsub_enum.c` and
`crsub_common.h`; rebuilding gives a binary with the recorded hash under gcc 13.3). The program needs x86-64 with BMI2 (`pext`/`pdep`).

Determinism: the stabiliser subgroups are generated from pseudo-random products, but every
random choice is derived from a fixed seed and the path from the root (`mix64` of the parent
seed and the orbit index, see `expand` and `dfs` in `search_n8/crsub_enum.c`). There is no external
randomness: a rerun produces the same tree, the same certificate (`cert8.sha256`) and the same
counts. Correctness does not depend on these choices (a smaller subgroup only gives more nodes).

Regenerating the certificate and the list of nodes of dimension 3 (in `search_n8/`):
`./crsub_enum 8 --good good8.bin --list3 --cert cert8.txt 2> list3_basis.txt`.

## 4. Expected numbers (`expected_counts.json`)

| Dimension | First program | Second program |
|---|---|---|
| 1 | 1 | (certificate) |
| 2 | 3 | (certificate) |
| 3 | 486 | 486 |
| 4 | 1,466,318 | 8,604,140 |
| 5 | 244,686,904 (244,668,896 fail P2; 18,008 fail P3) | 1,381,518,812 (1,381,416,881 fail P2; 101,931 fail P3) |
| 6 | 0 | 0 |

The numbers of the first program at dimensions 4 and 5 depend on the stabiliser subgroups it
finds (deterministically, see above); those of the second program depend only on the 486 nodes.

## 5. Dimension 6 (validation of the method)

In `search_n8/`, after building `crsub_enum`, `crsub_brute2` and `crsub_brute3`
(`gcc -O3 -march=native -o X X.c -lm`):

```
./crsub_enum 6 --out sol6.txt                  # 27 leaves, 9 distinct subspaces, 2 classes
./crsub_enum 6 --nosym --out sol6_nosym.txt    # 12,288 subspaces containing J
python3 crsub_invariants.py 6 sol6.txt --perm  # invariants of the two classes
./crsub_brute2 6 | sort -u | wc -l             # 12,288: enumeration without radical points (21 min)
./crsub_brute3 6 | cut -d' ' -f2- | sort -u | wc -l   # 12,288: constant rank only (5 min)
```
Logs: `logs/crsub_n6_sym.log`, `logs/crsub_brute2_n6.log`, `logs/crsub_brute3_n6.log`.

## 6. Small computations

`bash check_all.sh` (about a minute) re-checks the rank counts, the MacWilliams table, the SAT
model for n = 4 and n = 6, the Desarguesian instance for n = 8 and the examples for x^5.
It does **not** run the n = 8 search of Sections 2-3 above.

## 7. Environment of the reported runs

Ubuntu 24.04, two-core virtual machine, Intel Xeon 2.1 GHz; gcc 13.3 (`-O3 -march=native`);
Python 3.11, NumPy 2.4. First program: two parallel processes of 2.2 and 2.1 hours. Second
program: two parallel processes of 6.8 and 6.7 hours, partly on a machine shared with the second
run of the first program (a node that took 133 s then takes 67 s on an idle core).

## 8. Layout and logs

The files were sorted into folders on 2026-10-03, after all reported runs. The logs in `logs/`
and `second_impl/` were written with the earlier flat layout, so the command lines recorded in
them have no folder prefixes (for example `./crsub_enum`, `logs/...`). The two source files of
the first program were moved without change (same SHA-256 as in
`logs/crsub_n8_final_provenance.txt`). See `logs/NOTES.md` for the scripts that received path
lines.
