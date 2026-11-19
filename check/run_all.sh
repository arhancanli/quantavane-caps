#!/bin/sh
# run_all.sh: checks every cap of this repository with both checkers, in both layouts (caps/: the matrix layout of
# Edel's pages; points/: one point per line); makes sure that both checkers refuse the red control of every cap
# (check/red/, rewritten here by make_red.py: one point moved onto a line through two others) because of a collinear
# triple; computes the rank of every matrix and the sphere-packing bound (params.py); and recomputes the weight
# distributions written in caps/pg5_7_442.txt and caps/pg4_16_630.txt.
# usage: sh check/run_all.sh        (from the repository root; PYTHON = a python3 with numpy, CC = a C compiler)
# Writes check/results.txt. Exit 0 only if every cap passed both checkers in both layouts, both checkers refused
# every red control, every matrix has full rank, and both weight distributions matched.
set -u
PY=${PYTHON:-python3}
${CC:-cc} -O2 -o check/check_q_b check/check_q_b.c || exit 2
$PY -I check/make_red.py > /dev/null || exit 2
sha() { $PY -I -c 'import hashlib, sys; print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())' "$1"; }
bad=0
{
echo "run_all.sh $(date -u +%Y-%m-%dT%H:%M:%SZ); $($PY --version 2>&1); numpy $($PY -I -c 'import numpy; print(numpy.__version__)'); $(${CC:-cc} --version 2>&1 | head -1)"
while read name q model size; do
  echo
  echo "=== $name: $size points, q = $q, field model $model"
  for f in caps/$name.txt points/$name.txt; do
    echo "--- $f (sha256 $(sha $f))"
    $PY -I check/check_q_a.py $f $q $model $size; a=$?
    check/check_q_b $f $q $model $size; b=$?
    echo "exit codes: check_q_a.py $a, check_q_b $b"
    [ $a = 0 ] && [ $b = 0 ] || bad=1
  done
  r=check/red/$name.txt
  echo "--- red control $r (sha256 $(sha $r)): both checkers must refuse it"
  oa=$($PY -I check/check_q_a.py $r $q $model $size 2>&1); a=$?; echo "$oa"
  ob=$(check/check_q_b $r $q $model $size 2>&1); b=$?; echo "$ob"
  ra=0; rb=0
  case "$oa" in *"cap=False ("*"collinear"*) [ $a = 1 ] && ra=1;; esac
  case "$ob" in *"collinear_triple=1"*) [ $b = 1 ] && rb=1;; esac
  if [ $ra = 1 ] && [ $rb = 1 ]; then echo "red control: both checkers refused it (collinear triple), as required"
  else echo "red control: ERROR, a checker did not refuse it for a collinear triple"; bad=1; fi
done <<EOF
pg8_3_542 3 prime 542
pg10_3_2759 3 prime 2759
pg12_3_13428 3 prime 13428
pg9_5_17151 5 prime 17151
pg5_7_442 7 prime 442
pg8_7_21928 7 prime 21928
pg4_16_630 16 poly:11001 630
EOF
echo
echo "=== code parameters (check/params.py)"
$PY -I check/params.py caps/pg8_3_542.txt 3 prime || bad=1
$PY -I check/params.py caps/pg10_3_2759.txt 3 prime || bad=1
$PY -I check/params.py caps/pg12_3_13428.txt 3 prime || bad=1
$PY -I check/params.py caps/pg9_5_17151.txt 5 prime || bad=1
$PY -I check/params.py caps/pg5_7_442.txt 7 prime || bad=1
$PY -I check/params.py caps/pg8_7_21928.txt 7 prime || bad=1
$PY -I check/params.py caps/pg4_16_630.txt 16 poly:11001 || bad=1
echo
echo "=== weight distributions (check/weights.py)"
$PY -I check/weights.py caps/pg5_7_442.txt 7 prime || bad=1
$PY -I check/weights.py caps/pg4_16_630.txt 16 poly:11001 || bad=1
echo
if [ $bad = 0 ]; then echo "ALL PASSED"; else echo "SOMETHING FAILED"; fi
} > check/results.txt 2>&1
exit $bad
