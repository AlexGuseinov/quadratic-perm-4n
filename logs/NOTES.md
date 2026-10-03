# Notes on the logs

## constrank8_trinomial.log (first run, 2026-09-27)

The line `best bad=4` in this log is an artefact of early exit: `count_bad`
stops counting once the count exceeds `maxbad` (3 in this run), and the first
version of the program still reported that truncated value as a "best".
What the run does establish: among all 213,021,900 DO trinomials
x^d1 + c2 x^d2 + c3 x^d3 over GF(2^8), none has 3 or fewer components of
rank != 6 (histogram bins 0..3 are all zero). The true minimum was not
computed. Example check: x^3 + x^5 + x^9 actually has 127 bad components
(see trinomial_hyperplane.log). The reporting bug is fixed in constrank8.c.

## constrank8_f2coef.log

Exact: the counts in bins 0..3 are exact, so the minimum 3 (12288 functions)
is correct. For the example checked, the three bad components are exactly
b in F_4^* (f2coef_bad_components.log).

## constrank8_binomial.log

Exact (run with maxbad = 255): minimum 90.

## Folder reorganisation (2026-10-03)

All logs in this folder were written before the files were sorted into folders; command lines
and file names inside the logs refer to the earlier flat layout (everything in the repository
root). Nothing in `logs/` was rewritten.

- First program: `crsub_enum.c` and `crsub_common.h` were moved to `search_n8/` unchanged
  (same SHA-256 as in `crsub_n8_final_provenance.txt` and as the copies in
  `crsub_n8_final_src/`); the rebuilt binary has the recorded hash.
- Python scripts that import modules from another folder received two lines that add that
  folder to the import path, and scripts that write to `logs/` or call `bin/kissat` had those
  paths adjusted. No other line was changed. For this reason the hashes of `sym_search.py`,
  `sym_runner.py`, `sym_twostage.py`, `sym_case32.py`, `sym_case44.py` and `sym_p2.py` printed
  at the top of `sym_summary_n8.txt` differ from the current files; the table below the hashes
  is reproduced identically by `earlier_approaches/symmetry/sym_summary.py`.
- `du4perm_sat.py`: the hash in `sym_summary_n8.txt` (`f297b4d5...`) is that of the version
  used for the symmetric-case runs in the night of 2026-09-28; the file was edited later that
  morning (current hash `5934ff68...`), before the reorganisation. The earlier version was not
  kept.
