# gfq.py: small finite fields GF(q), q = p^k <= 32, as addition / multiplication tables on the symbols 0..q-1, and the
# encodings used to read Edel's generator matrices. Field model strings: "prime" (q prime, integers mod q),
# "poly:<c0c1...ck>" (symbol d = sum c_i p^i <-> polynomial sum c_i x^i modulo the monic polynomial with
# coefficients c0..ck, low to high) and "exp:<c0...ck>" (symbol 0 = 0, symbol d >= 1 = a^(d-1), a a root of that
# primitive polynomial). Used by qlns.py and check_q_a.py (check_q_b.c builds its own tables in C).
import itertools
import numpy as np


def factor(q):
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        k, r = 0, q
        while r % p == 0:
            r //= p; k += 1
        if r == 1 and k:
            return p, k
    raise ValueError(q)


def poly_tables(p, k, f):
    """symbol <-> polynomial over GF(p) of degree < k (digits base p, low first), arithmetic modulo monic f"""
    q = p ** k
    vec = [[(d // p ** i) % p for i in range(k)] for d in range(q)]
    def enc(v): return sum(c * p ** i for i, c in enumerate(v))
    def mul(a, b):
        r = [0] * (2 * k)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
        for t in range(2 * k - 1, k - 1, -1):
            c = r[t]
            if c:
                for j in range(k + 1):
                    r[t - k + j] = (r[t - k + j] - c * f[j]) % p
        return r[:k]
    ADD = np.array([[enc([(x + y) % p for x, y in zip(vec[a], vec[b])]) for b in range(q)] for a in range(q)], dtype=np.int16)
    MUL = np.array([[enc(mul(vec[a], vec[b])) for b in range(q)] for a in range(q)], dtype=np.int16)
    return ADD, MUL


def field(q, model):
    p, k = factor(q)
    if model == 'prime':
        assert k == 1
        a = np.arange(q)
        ADD = ((a[:, None] + a[None, :]) % q).astype(np.int16); MUL = ((a[:, None] * a[None, :]) % q).astype(np.int16)
    else:
        kind, coeffs = model.split(':')
        f = [int(c) for c in coeffs]; assert len(f) == k + 1 and f[-1] == 1
        ADD, MUL = poly_tables(p, k, f)
        if kind == 'exp':   # relabel: symbol d >= 1 means alpha^(d-1), alpha = the class of x
            x = p if k > 1 else None
            powers = [1]
            for _ in range(q - 2):
                powers.append(int(MUL[powers[-1], x]))
            assert len(set(powers)) == q - 1, 'polynomial is not primitive'
            lab = [0] + powers                      # symbol -> polynomial-basis value
            inv = {v: s for s, v in enumerate(lab)}
            ADD = np.array([[inv[int(ADD[lab[a], lab[b]])] for b in range(q)] for a in range(q)], dtype=np.int16)
            MUL = np.array([[inv[int(MUL[lab[a], lab[b]])] for b in range(q)] for a in range(q)], dtype=np.int16)
    # field axioms (cheap sanity): unique inverses, distributivity on all triples
    INV = np.zeros(q, dtype=np.int16)
    for a in range(1, q):
        hits = [b for b in range(1, q) if MUL[a, b] == 1]
        assert len(hits) == 1, f'{model}: no inverse for {a}'
        INV[a] = hits[0]
    for a, b, c in itertools.product(range(q), repeat=3):
        assert MUL[a, ADD[b, c]] == ADD[MUL[a, b], MUL[a, c]]
    return ADD, MUL, INV


def candidates(q):
    p, k = factor(q)
    if k == 1:
        return ['prime']
    out = []
    for tail in itertools.product(range(p), repeat=k):
        f = list(tail) + [1]
        if f[0] == 0:
            continue
        model = 'poly:' + ''.join(map(str, f))
        try:
            field(q, model)
        except AssertionError:
            continue                               # reducible
        out.append(model)
        try:
            field(q, 'exp:' + ''.join(map(str, f))); out.append('exp:' + ''.join(map(str, f)))
        except AssertionError:
            pass                                   # not primitive
    return out


SYMBOLS = '0123456789ABCDEFGHIJKLMNOPQRSTUV'


def read_matrix(path):
    """rows of an Edel generator-matrix page (or a plain matrix file): list of strings of symbols"""
    import re
    t = open(path, encoding='latin-1').read()
    rows = re.findall(r'^\s*([0-9A-Va-v]{3,})\s*(?:<br>|<BR>)?\s*$', t, re.M)
    n = len(rows[0]); assert all(len(r) == n for r in rows)
    return [r.upper() for r in rows]


def columns(rows):
    return [tuple(SYMBOLS.index(r[j]) for r in rows) for j in range(len(rows[0]))]
