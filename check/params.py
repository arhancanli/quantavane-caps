# params.py: the code parameters of a cap's matrix. The rank over GF(q) of the (k + 1) x m matrix (Gaussian
# elimination), so the code it checks has dimension m - rank; and the sphere-packing bound: if
# 1 + m (q - 1) + C(m, 2) (q - 1)^2 > q^rank, no code with this check matrix has minimum distance 5 or more. With the
# checkers' result (no three columns dependent, so minimum distance at least 4) this gives an [m, m - rank, 4]_q code.
# usage: python3 -I check/params.py file q model           (model as for check_q_a.py)
import sys, os
from math import comb
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gfq

path, q, model = sys.argv[1], int(sys.argv[2]), sys.argv[3]
ADD, MUL, INV = gfq.field(q, model)
rows = [[gfq.SYMBOLS.index(c) for c in l.strip().upper()] for l in open(path) if l.strip() and not l.startswith('#')]
k1, m = len(rows), len(rows[0])
M = np.array(rows, dtype=np.int64)
rank, col = 0, 0
while rank < k1 and col < m:
    nz = [i for i in range(rank, k1) if M[i, col]]
    if not nz:
        col += 1
        continue
    M[[rank, nz[0]]] = M[[nz[0], rank]]
    M[rank] = MUL[INV[M[rank, col]], M[rank]]
    for i in range(k1):
        if i != rank and M[i, col]:
            neg = next(t for t in range(q) if ADD[M[i, col], t] == 0)    # the additive inverse of M[i, col]
            M[i] = ADD[M[i], MUL[neg, M[rank]]]
    rank += 1
    col += 1
sphere = 1 + m * (q - 1) + comb(m, 2) * (q - 1) ** 2
print(f'P: file={path} q={q} matrix {k1} x {m}, rank {rank}; code length {m}, dimension {m - rank}; '
      f'sphere-packing: 1 + m(q-1) + C(m,2)(q-1)^2 = {sphere} > q^{rank} = {q ** rank}: {sphere > q ** rank}')
ok = rank == k1 and sphere > q ** rank
print(f'P: RESULT {"OK: with no three columns dependent, an [" + str(m) + "," + str(m - rank) + ",4]_" + str(q) + " code" if ok else "FAIL"}')
sys.exit(0 if ok else 1)
