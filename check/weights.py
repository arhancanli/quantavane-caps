# weights.py: the weight distribution A'_w of the code spanned by the rows of a cap's matrix (the dual of the code
# that the matrix checks), as Edel's pages give it for small caps. Every nonzero codeword is a nonzero multiple of a
# vector h, and its weight is the number of points P of the cap with h . P != 0, so each hyperplane h (scaled so that
# its first nonzero coordinate is 1) contributes q - 1 codewords of that weight; A'_0 = 1. It also prints the two
# identities sum_w A'_w = q^(k+1) and sum_w w A'_w = m (q - 1) q^k, and compares with the "# A'_w = c" lines
# of the file if it has them. It enumerates all (q^(k+1) - 1)/(q - 1) hyperplanes, so it is for small spaces only.
# usage: python3 -I check/weights.py file q model           (model as for check_q_a.py)
import sys, os, re, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gfq

path, q, model = sys.argv[1], int(sys.argv[2]), sys.argv[3]
ADD, MUL, INV = gfq.field(q, model)
text = open(path).read().splitlines()
rows = [l.strip().upper() for l in text if l.strip() and not l.startswith('#')]
P = np.array(gfq.columns(rows), dtype=np.int64)          # m x (k + 1)
m, n = P.shape
hyper = np.array([v for v in itertools.product(range(q), repeat=n) if any(v) and [x for x in v if x][0] == 1],
                 dtype=np.int64)
assert len(hyper) == (q ** n - 1) // (q - 1)
dist = {0: 1}
for s in range(0, len(hyper), 4096):
    Hc = hyper[s:s + 4096]
    acc = np.zeros((len(Hc), m), dtype=np.int64)
    for i in range(n):
        acc = ADD[acc, MUL[Hc[:, i][:, None], P[None, :, i]]]
    for w, c in zip(*np.unique((acc != 0).sum(axis=1), return_counts=True)):
        dist[int(w)] = dist.get(int(w), 0) + int(c) * (q - 1)
dist = dict(sorted(dist.items()))
total = sum(dist.values())
first = sum(w * c for w, c in dist.items())
print(f"W: file={path} PG({n - 1},{q}) points={m} hyperplanes={len(hyper)}")
print("W: " + ", ".join(f"A'_{w} = {c}" for w, c in dist.items()))
print(f"W: sum A'_w = {total} (q^(k+1) = {q ** n}: {total == q ** n}); sum w A'_w = {first} "
      f"(m (q-1) q^k = {m * (q - 1) * q ** (n - 1)}: {first == m * (q - 1) * q ** (n - 1)})")
print(f"W: largest weight {max(dist)}; weight of the last row {int((P[:, -1] != 0).sum())}")
written = {int(w): int(c) for w, c in re.findall(r"A'_(\d+) = (\d+)", '\n'.join(l for l in text if l.startswith('#')))}
if written:
    print(f"W: the distribution written in the file {'matches' if written == dist else 'DOES NOT MATCH'}")
ok = total == q ** n and first == m * (q - 1) * q ** (n - 1) and (not written or written == dist)
print('W: RESULT', 'OK' if ok else 'FAIL')
sys.exit(0 if ok else 1)
