# check_q_a.py (checker A for caps in PG(n-1,q), Python/numpy; written separately from check_q_b.c).
# Reads a point file (one point per line, n symbols 0-9A-V; '#' lines are comments), a matrix file in the layout of
# Edel's pages (one row per coordinate, one column per point; recognised by rows longer than 32 symbols) or an Edel
# matrix page (.html, columns are points), normalises every point (first nonzero coordinate 1), and checks: symbols < q, no repeated
# projective point, the stated size, and that no three points are collinear: for every pair u, v the q - 1 other
# points of their line, u + l v (l != 0), are looked up in a table of the cap's normalised points.
# --red: red control, adds the point <u + v> of the first two points (must FAIL).
# usage: python -I check_q_a.py file q model expected_size [--red]      (model: prime | poly:.. | exp:.., see gfq.py)
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gfq

path, q, model, expected = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
ADD, MUL, INV = gfq.field(q, model)
if path.endswith('.html'):
    pts = gfq.columns(gfq.read_matrix(path))
else:
    lines = [l.strip().upper() for l in open(path) if l.strip() and not l.startswith('#')]
    if len(lines[0]) > 32:   # matrix layout, as on Edel's pages: one row per coordinate, one column per point
        if len(set(map(len, lines))) != 1:
            print(f'A: file={path} format=False (matrix rows of different lengths)')
            print('A: RESULT FAIL')
            sys.exit(1)
        pts = gfq.columns(lines)
    else:
        pts = [tuple(gfq.SYMBOLS.index(c) for c in l) for l in lines]
n = len(pts[0])
ok_fmt = all(len(p) == n and max(p) < q for p in pts) and all(any(p) for p in pts)
P = np.array(pts, dtype=np.int64)

def normalise(V):
    lead = np.zeros(len(V), dtype=np.int64)
    for j in range(V.shape[1] - 1, -1, -1):
        lead = np.where(V[:, j] != 0, V[:, j], lead)
    return MUL[INV[lead][:, None], V].astype(np.int64)

pw = q ** np.arange(n, dtype=np.int64)
def is_cap(P):
    N = normalise(P); idx = N @ pw
    if len(np.unique(idx)) != len(idx):
        return False, 'repeated projective point'
    member = np.zeros(q ** n, dtype=bool); member[idx] = True
    for i in range(len(N) - 1):
        V = N[i + 1:]
        for l in range(1, q):
            W = ADD[N[i][None, :], MUL[l][V]]
            hit = np.nonzero(member[normalise(W) @ pw])[0]
            if len(hit):
                return False, f'points {i}, {i + 1 + int(hit[0])} and a third are collinear'
    return True, ''

ok, why = is_cap(P)
print(f'A: file={path} PG({n - 1},{q}) field={model} points={len(P)} format={ok_fmt} size==expected({expected})={len(P) == expected} '
      f'cap={ok}' + (f' ({why})' if why else ''))
good = ok_fmt and ok and len(P) == expected
if '--red' in sys.argv:
    N = normalise(P[:2]); w = ADD[N[0], N[1]]
    ok_red, why = is_cap(np.vstack([P, w[None, :]]))
    print(f'A: red control (added the point <u+v> of the first two points): cap={ok_red} -> '
          f'{"FAIL as required" if not ok_red else "ERROR: control passed"}')
    good = good and not ok_red
print('A: RESULT', 'PASS' if good else 'FAIL')
sys.exit(0 if good else 1)
