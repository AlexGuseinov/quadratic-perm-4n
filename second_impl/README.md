# Second implementation of the n = 8 search

Written in Python/numpy by a separate agent from the mathematical specification of the search
only (it did not see `crsub_enum.c`). Data files are read from `$N8_DATA` (default:
`second_impl/data/`).

| File | Purpose |
|---|---|
| `task1_pfaffian.py` | Rank-6 table of all 2^28 forms by Pfaffians (bit-sliced); writes `data/rank6_pf.bin`, compares with `data/good8.bin` (the C table). Log: `task1.log` (149,920,960 forms; identical) |
| `verify_cert.py` | Checks the certificate of the symmetry reductions at levels 1-2 (`data/cert8.txt`, from `crsub_enum 8 --list3 --cert`) against candidate sets computed from the definitions, and the 486 roots (`data/list3_basis.txt`). Log: `verify_cert.log` (all checks OK) |
| `search.py` | Searches the subtree of a root (dimension-3 node) without symmetry; batch mode over all roots (`--batch data/list3_basis.txt --out RESULTS --start I --step M`, resumable) |
| `search_v1.py` | Earlier, slower version (produced `root0.log`, `root5.log`, `root17.log`) |
| `indep_lib.py`, `selftest.py`, `rootprobe.py` | Primitives and self-tests |
| `run_all.sh` | One command for the whole chain: Pfaffian table, certificate check, search of the 486 subtrees, summary (see `../ARTIFACT.md`) |
| `summarize_batch.py` | Checks that every root was processed with no solution and compares the node counts with `../expected_counts.json` |
| `batch_results_*.txt` | Results of the complete batch run over the 486 roots |

The certificate is shipped compressed (`data/cert8.txt.xz`, 1,468,672 lines; SHA-256 of the
uncompressed file in `data/cert8.sha256`); `verify_cert.py` reads either form. The rank tables
(32 MB each) are not shipped: `data/good8.bin` is rebuilt by the command below and
`data/rank6_pf.bin` by `task1_pfaffian.py`.

Regenerating the data: `search_n8/crsub_enum 8 --good second_impl/data/good8.bin --list3 --cert
second_impl/data/cert8.txt 2> second_impl/data/list3_basis.txt` (from the repository root; the
listing and the certificate are deterministic), then `python3 task1_pfaffian.py`.
