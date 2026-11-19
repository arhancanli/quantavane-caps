/* check_q_b.c (checker B for caps in PG(n-1,q), C; written separately from check_q_a.py).
   Field: built here from the model string: "prime" (integers mod q), "poly:c0..ck" (symbol = sum c_i p^i, arithmetic
   of polynomials over GF(p) modulo the monic polynomial c0 + c1 x + ... + ck x^k) or "exp:c0..ck" (symbol 0 = 0,
   symbol d = a^(d-1), a a root of that primitive polynomial).
   Reads a point file (one point per line, symbols 0-9A-V, '#' comments) or a matrix file in the layout of Edel's
   pages (one row per coordinate, one column per point; recognised by rows longer than 32 symbols), scales every
   point so that its first nonzero coordinate is 1, packs it into a 64-bit code (ceil(log2 q) bits per coordinate), sorts, and checks: no duplicate,
   the stated size, and for every pair a, b and every nonzero l that the point l*a + b (scaled) is not in the cap
   (binary search). With -r: red control, appends the scaled point a + b of the first two points (must FAIL).
   usage: check_q_b file q model expected_size [-r]          build: cc -O2 -o check_q_b check_q_b.c            */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned long long u64;
static int q, p, k, n, BITS, ADD[32][32], MUL[32][32], INV[32];

static void build_field(const char *model) {
  int t = q; p = 0; k = 0;
  for (int c = 2; c <= 31 && !p; c++) if (t % c == 0) { p = c; while (t % c == 0) { t /= c; k++; } }
  if (!strcmp(model, "prime")) {
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) { ADD[a][b] = (a + b) % q; MUL[a][b] = (a * b) % q; }
  } else {
    const char *c = strchr(model, ':') + 1; int f[8]; for (int i = 0; i <= k; i++) f[i] = c[i] - '0';
    int pa[32][32], pm[32][32];  /* polynomial-basis tables */
    for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) {
      int x[8], y[8], s = 0, w = 1, r[16] = {0};
      for (int i = 0, u = a, v = b; i < k; i++, u /= p, v /= p) { x[i] = u % p; y[i] = v % p; }
      for (int i = 0; i < k; i++) { s += ((x[i] + y[i]) % p) * w; w *= p; }
      pa[a][b] = s;
      for (int i = 0; i < k; i++) for (int j = 0; j < k; j++) r[i + j] = (r[i + j] + x[i] * y[j]) % p;
      for (int d = 2 * k - 2; d >= k; d--) if (r[d]) { int cc = r[d]; for (int j = 0; j <= k; j++) r[d - k + j] = ((r[d - k + j] - cc * f[j]) % p + p) % p; }
      s = 0; w = 1; for (int i = 0; i < k; i++) { s += r[i] * w; w *= p; }
      pm[a][b] = s;
    }
    if (!strncmp(model, "poly", 4)) {
      memcpy(ADD, pa, sizeof pa); memcpy(MUL, pm, sizeof pm);
    } else {  /* exp: symbol d >= 1 is x^(d-1) */
      int lab[32], inv[32]; lab[0] = 0; lab[1] = 1;
      for (int d = 2; d < q; d++) lab[d] = pm[lab[d - 1]][p];
      for (int d = 0; d < q; d++) inv[lab[d]] = d;
      for (int a = 0; a < q; a++) for (int b = 0; b < q; b++) { ADD[a][b] = inv[pa[lab[a]][lab[b]]]; MUL[a][b] = inv[pm[lab[a]][lab[b]]]; }
    }
  }
  for (int a = 1; a < q; a++) { INV[a] = -1; for (int b = 1; b < q; b++) if (MUL[a][b] == 1) INV[a] = b; if (INV[a] < 0) { fprintf(stderr, "not a field\n"); exit(2); } }
}
static u64 pack(const int *v) { u64 c = 0; for (int i = 0; i < n; i++) c |= (u64)v[i] << (BITS * i); return c; }
static void unpack(u64 c, int *v) { for (int i = 0; i < n; i++) v[i] = (int)((c >> (BITS * i)) & ((1u << BITS) - 1)); }
static void scale(int *v) { int lead = 0; for (int i = 0; i < n; i++) if (v[i]) { lead = v[i]; break; } if (!lead) return; int s = INV[lead]; for (int i = 0; i < n; i++) v[i] = MUL[s][v[i]]; }
static int cmp(const void *a, const void *b) { u64 x = *(const u64 *)a, y = *(const u64 *)b; return x < y ? -1 : x > y; }
static int has(const u64 *s, long m, u64 c) { long lo = 0, hi = m - 1; while (lo <= hi) { long mid = (lo + hi) / 2; if (s[mid] == c) return 1; if (s[mid] < c) lo = mid + 1; else hi = mid - 1; } return 0; }
static long collinear(u64 *pts, long m) {   /* number of pairs whose line holds a third cap point (stops at 1) */
  u64 *s = malloc(m * sizeof(u64)); memcpy(s, pts, m * sizeof(u64)); qsort(s, m, sizeof(u64), cmp);
  int a[32], b[32], w[32];
  for (long i = 0; i < m; i++) { unpack(pts[i], a);
    for (long j = i + 1; j < m; j++) { unpack(pts[j], b);
      for (int l = 1; l < q; l++) { for (int t = 0; t < n; t++) w[t] = ADD[MUL[l][a[t]]][b[t]]; scale(w); if (has(s, m, pack(w))) { free(s); return 1; } } } }
  free(s); return 0;
}
int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: check_q_b file q model expected_size [-r]\n"); return 2; }
  q = atoi(argv[2]); build_field(argv[3]); BITS = 1; while ((1 << BITS) < q) BITS++; long expected = atol(argv[4]); int red = argc > 5 && !strcmp(argv[5], "-r");
  FILE *f = fopen(argv[1], "r"); if (!f) { perror("open"); return 2; }
  const char *SYM = "0123456789ABCDEFGHIJKLMNOPQRSTUV";
  /* every non-comment line is kept; if the first one is longer than 32 symbols the file is a matrix in the layout of
     Edel's pages (one row per coordinate, one column per point), otherwise a point file (one point per line) */
  char *line = NULL; size_t lcap = 0; ssize_t got; char **rows = NULL; long *rlen = NULL, nr = 0, rcap = 0;
  while ((got = getline(&line, &lcap, f)) != -1) {
    if (line[0] == '#') continue;
    long len = 0; while (line[len] && line[len] != '\n' && line[len] != '\r') len++;
    if (!len) continue;
    if (nr == rcap) { rcap = rcap ? 2 * rcap : 64; rows = realloc(rows, rcap * sizeof(char *)); rlen = realloc(rlen, rcap * sizeof(long)); }
    rows[nr] = malloc(len + 1); memcpy(rows[nr], line, len); rows[nr][len] = 0; rlen[nr++] = len;
  }
  fclose(f); free(line);
  int matrix = nr > 0 && rlen[0] > 32;
  long m = 0, cap = 1024, npts = nr == 0 ? 0 : matrix ? rlen[0] : nr; u64 *pts = malloc(cap * sizeof(u64)); int fmt = 1;
  n = nr == 0 ? 0 : matrix ? (int)nr : (int)rlen[0];
  if (n < 1 || n > 32 || n * BITS > 64) { fmt = 0; npts = 0; }
  if (matrix) for (long i = 0; i < nr; i++) if (rlen[i] != rlen[0]) { fmt = 0; npts = 0; }
  for (long j = 0; j < npts; j++) {
    if (!matrix && rlen[j] != n) { fmt = 0; continue; }
    int v[32], nz = 0;
    for (int i = 0; i < n; i++) { char ch = matrix ? rows[i][j] : rows[j][i]; if (ch >= 'a' && ch <= 'v') ch -= 32; const char *pp = strchr(SYM, ch); if (!pp || (pp - SYM) >= q) { fmt = 0; v[i] = 0; } else v[i] = (int)(pp - SYM); nz |= v[i]; }
    if (!nz) fmt = 0;
    scale(v);
    if (m == cap) { cap *= 2; pts = realloc(pts, cap * sizeof(u64)); }
    pts[m++] = pack(v);
  }
  u64 *s = malloc(m * sizeof(u64)); memcpy(s, pts, m * sizeof(u64)); qsort(s, m, sizeof(u64), cmp);
  int dup = 0; for (long i = 1; i < m; i++) if (s[i] == s[i - 1]) dup = 1;
  long col = collinear(pts, m);
  int good = fmt && !dup && m == expected && !col;
  printf("B: PG(%d,%d) field=%s points=%ld format=%d duplicates=%d size==expected(%ld)=%d collinear_triple=%ld\n", n - 1, q, argv[3], m, fmt, dup, expected, m == expected, col);
  if (red && m >= 2) {
    int a[32], b[32], w[32]; unpack(pts[0], a); unpack(pts[1], b);
    for (int t = 0; t < n; t++) w[t] = ADD[a[t]][b[t]]; scale(w);
    pts = realloc(pts, (m + 1) * sizeof(u64)); pts[m] = pack(w);
    long c2 = collinear(pts, m + 1);
    printf("B: red control (point <a+b> of the first two points appended): collinear_triple=%ld -> %s\n", c2, c2 ? "FAIL as required" : "ERROR: control passed");
    good = good && c2;
  }
  printf("B: RESULT %s\n", good ? "PASS" : "FAIL");
  return good ? 0 : 1;
}
