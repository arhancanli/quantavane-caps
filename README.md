# Caps in PG(8,3), PG(10,3), PG(12,3), PG(9,5), PG(5,7), PG(8,7) and PG(4,16)

A cap in the projective space PG(k,q) is a set of points no three of which are collinear. This repository holds
seven caps found by Quantavane in November 2026. Each is larger than the cap listed for its space in Yves Edel's
tables of large caps (www.yvesedel.de/Matritzen/CAPs), and each was found from that cap, in its coordinates. Every
cap is given as a check matrix in the layout of Edel's pages and as a list of points, with two separately written
checking programs, their output, and a red control for each cap.

| Space | Size | Previous record (source) | Edel's cap that the search started from | File |
|---|---|---|---|---|
| PG(8,3) | 542 | 541 (Edel, yvesedel.de and Bierbrauer-Edel survey, Table 1 (2012)) | 533 of its 541 points kept, 9 new points | `caps/pg8_3_542.txt` |
| PG(10,3) | 2759 | 2744 (Edel, yvesedel.de and Bierbrauer-Edel survey, Table 1 (2012)) | 2718 of its 2744 points kept, 41 new points | `caps/pg10_3_2759.txt` |
| PG(12,3) | 13428 | 13312 (Edel, yvesedel.de (2010)) | 13143 of its 13312 points kept, 285 new points | `caps/pg12_3_13428.txt` |
| PG(9,5) | 17151 | 17124 (Edel, yvesedel.de and Bierbrauer-Edel survey, Table 1 (2012)) | 17073 of its 17124 points kept, 78 new points | `caps/pg9_5_17151.txt` |
| PG(5,7) | 442 | 434 (Edel, yvesedel.de and Bierbrauer-Edel survey, Table 1 (2012)) | 410 of its 434 points kept, 32 new points | `caps/pg5_7_442.txt` |
| PG(8,7) | 21928 | 21555 (Edel, yvesedel.de and Bierbrauer-Edel survey, Table 1 (2012)) | 20461 of its 21555 points kept, 1467 new points | `caps/pg8_7_21928.txt` |
| PG(4,16) | 630 | 629 (Edel, yvesedel.de (2001/2010)) | 616 of its 629 points kept, 14 new points | `caps/pg4_16_630.txt` |

Edel's index page of large caps, read on 19 November 2026, was last changed on 28 April 2010. On 12 November 2026
we read Edel's tables and the Bierbrauer–Edel survey for each of these spaces and searched the literature; we
found no larger cap in any of them. These searches were limited, so novelty is our reading of them, not an
established fact.

## As codes

The matrix of a cap of size m in PG(k,q) is a check matrix of a linear code of length m. `check/params.py` finds that
each matrix has rank k + 1, so the code has dimension m − k − 1. Its minimum distance is at least 4 because no three
columns are linearly dependent (the checkers below), and at most 4 by the sphere-packing bound (also checked by
`check/params.py`).

## Files

- `caps/<name>.txt`: each cap as a check matrix in the layout of Edel's pages (below).
- `points/<name>.txt`: the same points, one per line (coordinates x0, ..., xk), in the order of the matrix columns.
- `check/check_q_a.py` with `check/gfq.py`: checker A (Python with numpy).
- `check/check_q_b.c`: checker B (C).
- `check/make_red.py` and `check/red/<name>.txt`: the red controls.
- `check/params.py`: rank of each matrix and the sphere-packing bound.
- `check/weights.py`: the weight distributions given for the two smaller caps.
- `check/run_all.sh`: runs all of the above; `check/results.txt` is its output.
- `SHA256SUMS`: the SHA-256 of every other file.
- `LICENSE` (CC BY 4.0), `LICENSE-CODE` (MIT), `CITATION.cff`.

### The matrix files

As on Edel's pages, the points of a cap of size m in PG(k,q) are the columns of a (k + 1) × m matrix,
written as k + 1 rows of m symbols each, with no separators; the first row holds x0, the second x1, and so on.
Over GF(3), GF(5) and GF(7) the symbols are the integers 0, ..., q − 1, with arithmetic modulo q. Over GF(16) the
encoding is Edel's: the prime polynomial used to generate GF(16) is X^4 + X + 1, and the element
aX^3 + bX^2 + cX + d (a, b, c, d in {0, 1}) is written as the number 8a + 4b + 2c + d, in hexadecimal (0–9, A–F).
Lines that start with `#` are comments: what the cap is, the field, and the order of the columns.

The coordinates are those of Edel's matrix for the cap that the search started from, so the columns can be compared
with his directly. Every column is scaled so that its last nonzero coordinate is 1. The columns come in two groups:
first the points of Edel's cap that were kept, then the new points, each group in lexicographic order of
(x0, ..., xk); the comment lines give the column numbers of the two groups. On Edel's pages the last row of each
matrix is a word of maximal weight in the code spanned by the rows. Because the coordinates are kept, that is not
always so here: for the 442-cap the last row has weight 434, while the largest weight is 442 (one hyperplane,
x0 + x4 + 3x5 = 0, contains no point of the cap); for the 630-cap the last row has weight 630, the largest.

For the two smaller caps, in PG(5,7) and PG(4,16), the file ends with the weight distribution A'_w of the code
spanned by the rows, as Edel's pages give it for small caps. `check/weights.py` computes it from the matrix by going
through all hyperplanes.

## How the caps were found

By exact-refill large-neighbourhood search, started from the cap in Edel's tables for the same space. Each step
removes a few dozen points (ends of the secants through outside points that lie on few secants, plus random points),
collects every point whose secants all meet the removed set, and lets the CP-SAT solver (Google OR-Tools 9.15) find
the largest cap among them; a larger refill is a larger cap. For several of the caps, the points added by different
runs were pooled with Edel's cap and solved exactly once more. The search code is not part of this repository; the
caps do not depend on it, since they are checked directly.

## How to check

You need Python 3 with numpy and a C compiler.

```sh
sh check/run_all.sh        # both checkers on every cap in both layouts, the red controls, params and weights
                           # (PYTHON=python3 and CC=cc by default); writes check/results.txt
python3 -I check/check_q_a.py caps/pg8_3_542.txt 3 prime 542
cc -O2 -o check/check_q_b check/check_q_b.c
check/check_q_b caps/pg8_3_542.txt 3 prime 542
```

The arguments are the file, q, the field model (`prime`, or `poly:11001` for GF(16) = GF(2)[x]/(x^4 + x + 1) in
Edel's encoding) and the expected number of points. Each checker reads either layout (a file whose rows are longer
than 32 symbols is read as a matrix) and prints `RESULT PASS` or `RESULT FAIL`.

- Checker A (`check_q_a.py`, Python with numpy): scales every point so that its first nonzero coordinate is 1;
  checks that the symbols are field elements, that no point is zero, that no projective point is repeated and that
  the size is as stated; then, for every two points u and v, looks up the other q − 1 points u + λv (λ ≠ 0) of their
  line in a table of the cap's points.
- Checker B (`check_q_b.c`, C, with its own field tables): the same scaling and the same checks of the input, with
  the points as sorted 64-bit codes; then, for every two points a and b and every λ ≠ 0, a binary search for the point
  λa + b.
- Red controls (`check/make_red.py`): for a cap with m points, the column m//2 + 1 is replaced by the sum of the
  columns m//3 + 1 and 2m//3 + 1 (scaled), a point on the line through those two points. The damaged file has the
  same size and no repeated point. `run_all.sh` requires both checkers to refuse it because of a collinear triple.

Both checkers were written on 12 November 2026 for the search that found these caps, and used there to check them.
They are unchanged here except that each also reads the matrix layout. Both still have their built-in red control
(`--red` for A, `-r` for B: the point u + v of the first two points is added); `run_all.sh` uses the red-control
files instead.

## Results

`check/results.txt` is the output of `sh check/run_all.sh` on 19 November 2026 (Python 3.12.13, numpy 2.5.4, Apple
clang 21.0.0):

| Cap | Checker A: matrix, points | Checker B: matrix, points | Red control: A, B |
|---|---|---|---|
| 542 in PG(8,3) | PASS, PASS | PASS, PASS | refused, refused |
| 2759 in PG(10,3) | PASS, PASS | PASS, PASS | refused, refused |
| 13428 in PG(12,3) | PASS, PASS | PASS, PASS | refused, refused |
| 17151 in PG(9,5) | PASS, PASS | PASS, PASS | refused, refused |
| 442 in PG(5,7) | PASS, PASS | PASS, PASS | refused, refused |
| 21928 in PG(8,7) | PASS, PASS | PASS, PASS | refused, refused |
| 630 in PG(4,16) | PASS, PASS | PASS, PASS | refused, refused |

Every red control was refused by both checkers because of a collinear triple. Every matrix has rank k + 1, and the
sphere-packing bound excludes minimum distance 5 for all seven codes (`check/params.py`). The weight distributions
recomputed by `check/weights.py` match those written in `caps/pg5_7_442.txt` and `caps/pg4_16_630.txt`.

## Authorship and the role of AI

Author and maintainer: Arhan Canli (Quantavane, https://quantavane.com), who takes responsibility for the results.

The searches, the checking programs and the runs were written and carried out by AI agents (Claude, by Anthropic)
working for Quantavane under Arhan Canli's direction. Programs described as written separately share no code but
were written with the same AI system. No mathematician outside Quantavane has reviewed this repository; the caps can
be checked by anyone with the programs here or with their own.

## How to cite

If you use these caps, please cite them as below, with Yves Edel's tables of large caps, from whose caps they were
found.

```bibtex
@misc{quantavane-caps,
  author       = {Canli, Arhan},
  title        = {Caps in {PG}(8,3), {PG}(10,3), {PG}(12,3), {PG}(9,5), {PG}(5,7), {PG}(8,7) and {PG}(4,16)},
  year         = {2026},
  howpublished = {\url{https://github.com/arhancanli/quantavane-caps}},
  note         = {Quantavane; CC BY 4.0 for the data and text, MIT licence for the code}
}
```

[CITATION.cff](CITATION.cff) gives the same reference for GitHub's "Cite this repository" button.

## Licence

The caps, the red-control files, the output and the text are under CC BY 4.0 (`LICENSE`); the code is under the MIT
licence (`LICENSE-CODE`).

## References

- Y. Edel, tables of large caps with their generator matrices, www.yvesedel.de/Matritzen/CAPs/CAPMatIndex.html
  (formerly at www.mathi.uni-heidelberg.de/~yves/), last changed 28 April 2010.
- J. Bierbrauer and Y. Edel, Large caps in projective Galois spaces, in: J. De Beule and L. Storme (eds.), Current
  Research Topics in Galois Geometry, Nova Science Publishers, 2012, 87–104.
