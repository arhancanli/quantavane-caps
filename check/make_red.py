# make_red.py: red controls. For every cap, writes check/red/<name>.txt: the cap's matrix with one point moved onto
# a line through two others. With m points P_1 ... P_m (the columns of caps/<name>.txt), the column at position
# m // 2 + 1 is replaced by P_a + P_b (scaled so that its last nonzero coordinate is 1), where a = m // 3 + 1 and
# b = 2m // 3 + 1. The damaged file has the same size and no repeated point, but three collinear points, so both
# checkers must refuse it.
# usage: python3 -I check/make_red.py          (from the repository root)
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gfq

CAPS = [('pg8_3_542', 3, 'prime'), ('pg10_3_2759', 3, 'prime'), ('pg12_3_13428', 3, 'prime'),
        ('pg9_5_17151', 5, 'prime'), ('pg5_7_442', 7, 'prime'), ('pg8_7_21928', 7, 'prime'),
        ('pg4_16_630', 16, 'poly:11001')]
os.makedirs('check/red', exist_ok=True)
for name, q, model in CAPS:
    ADD, MUL, INV = gfq.field(q, model)
    lines = [l.rstrip('\r\n') for l in open(f'caps/{name}.txt')]
    rows = [l for l in lines if l and not l.startswith('#')]
    m = len(rows[0])
    cols = [[gfq.SYMBOLS.index(r[j]) for r in rows] for j in range(m)]

    def scaled(v):
        s = INV[[x for x in v if x][-1]]
        return tuple(int(MUL[s, x]) for x in v)
    a, b, r = m // 3, 2 * m // 3, m // 2
    w = scaled([int(ADD[x, y]) for x, y in zip(cols[a], cols[b])])
    assert w not in set(scaled(c) for c in cols)
    cols[r] = list(w)
    with open(f'check/red/{name}.txt', 'w') as f:
        f.write(f'# RED CONTROL (not a cap): caps/{name}.txt with column {r + 1} replaced by the sum of columns {a + 1} and '
                f'{b + 1}\n# (scaled), a point on the line through those two points. Both checkers must refuse this file.\n')
        f.write('\n'.join(''.join(gfq.SYMBOLS[c[i]] for c in cols) for i in range(len(rows))) + '\n')
    print(f'{name}: column {r + 1} replaced by column {a + 1} + column {b + 1} = {"".join(gfq.SYMBOLS[x] for x in w)}')
