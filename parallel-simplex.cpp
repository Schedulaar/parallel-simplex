#include <cstdio>
#include <cmath>
#include <string>
#include <bsp.hpp>
#include <random>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <chrono>
#include <limits>

using flt = double;
flt F_MAX = std::numeric_limits<flt>::infinity();
flt EPS = 1e-10;

enum STATUS {
  RUNNING, SUCCESS, UNBOUNDED
};

struct result {
  double z;
  long iters;
  STATUS status;
};

long M, N;
long nArg = -1;
int Gargc;
char **Gargv;

bool PRINT_TABLES = false;
bool PROFILING = false;
bool DEBUG = false;

double lastTime;
int NUM_STEPS = 9;
double times[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

void stepFinished(long step, long s, long t) {
  if (s == 0 && t == 0) {
    double newTime = bsp_time();
    times[step] += newTime - lastTime;
    lastTime = newTime;
  }
}

long min(long a, long b) { return a <= b ? a : b; }

void bsp_broadcast(flt *x, long n, long src, long s0, long s, long stride, long p0, long phase) {
  /* Broadcast the vector x of length n from processor src to
     processors s0+t*stride, 0 <= t < p0. Here n, p0 >= 1.
     The vector x must have been registered previously.
     Processors are numbered in 1D fashion.

     phase = phase of two-phase broadcast (0 or 1)
     Only one phase is performed, without synchronization.
  */

  if (n < 1 || p0 < 1)
    return;

  long b = (n % p0 == 0 ? n / p0 : n / p0 + 1); // broadcast block size

  if ((phase == 0 && s == src) ||
      (phase == 1 && s0 <= s && s < s0 + p0 * stride && (s - s0) % stride == 0)) {
    /* Participate */

    for (long t = 0; t < p0; t++) {
      long dest = s0 + t * stride;
      long t1 = (phase == 0 ? t : (s - s0) / stride);
      // in phase 1: s= s0+t1*stride

      long nbytes = min(b, n - t1 * b);
      if (nbytes > 0 && dest != src)
        bsp_put(dest, &x[t1 * b], x, t1 * b * sizeof(flt), min(b, n - t1 * b) * sizeof(flt));
    }
  }
}

flt **matallocd(size_t m, size_t n) {
  /* This function allocates an m x n matrix of flts */
  size_t i;
  flt *pd, **ppd;
  ppd = new flt *[m];
  if (ppd == NULL)
    bsp_abort("matallocd: not enough memory");
  pd = new flt[m * n];
  if (pd == NULL)
    bsp_abort("matallocd: not enough memory");
  ppd[0] = pd;
  for (i = 1; i < m; i++)
    ppd[i] = ppd[i - 1] + n;

  return ppd;
}

void swap(long &a, long &b) {
  long tmp = a;
  a = b;
  b = tmp;
}

long nloc(long p, long s, long n) {
  /* Compute number of local components of processor s for vector
     of length n distributed cyclically over p processors.
     Also useful for block distribution with
          ceil(n/p)  components if s < (n mod p),
          floor(n/p) components otherwise.
  */
  return (n + p - s - 1) / p;
}

void matfreed(flt **ppd) {
  /* This function frees a matrix of flts */
  if (ppd != NULL) {
    if (ppd[0] != NULL)
      free(ppd[0]);
    free(ppd);
  }
}

void print_tableau(long s, long t, long m, long n, flt **A, flt *c, flt *b, flt v,
                   long *Basis, long *NonBasis) {
  if (s == 0 && t == 0) {
    printf("--------------------\n");
    printf("z  = %+.2f", v);
  }
  bsp_sync();
  for (int j = 0; j < n; j++) {
    if (j % N == t && s == 0)
      printf(" %+.2f*x%li", c[j / N], NonBasis[j]);
    bsp_sync();
  }
  if (s == 0 && t == 0)
    printf("\n");
  bsp_sync();
  for (int i = 0; i < m; i++) {
    if (i % M == s && t == 0)
      printf("x%li = %+.2f", Basis[i], b[i / M]);
    bsp_sync();
    for (int j = 0; j < n; j++) {
      if (i % M == s && j % N == t)
        printf(" %+.2f*x%li", -A[i / M][j / N], NonBasis[j]);
      bsp_sync();
    }
    if (s == 0 && t == 0)
      printf("\n");
    bsp_sync();
  }
}

struct IndexValuePair {
  long index;
  flt value;
};

void findColPhase1(long s, long t, long locCols, flt *c, IndexValuePair *cLocalMaxima) {
  if (s == 0) {
    long localMaxInd = -1;
    flt localMaxVal = 0;
    for (long j = 0; j < locCols; j++) {
      if (c[j] > localMaxVal) {
        localMaxInd = j;
        localMaxVal = c[j];
      }
    }

    // Broadcast local maximum, if there is a possible index
    cLocalMaxima[t] = {.index = t + N * localMaxInd, .value = localMaxVal};
    for (long k = 0; k < N; k++) {
      if (k == t) continue;
      cLocalMaxima[k] = {.index=-1, .value=0};
      if (localMaxVal > EPS)
        bsp_put(0 * N + k, &cLocalMaxima[t], cLocalMaxima, t * sizeof(IndexValuePair), sizeof(IndexValuePair));
    }
  }
}

void findColPhase2(long s, long t, IndexValuePair *cLocalMaxima, long &e, STATUS &status) {
  if (s == 0) {
    long globalMaximumProc = 0;
    for (long k = 1; k < N; k++) {
      if (cLocalMaxima[globalMaximumProc].value < cLocalMaxima[k].value)
        globalMaximumProc = k;
    }
    e = cLocalMaxima[globalMaximumProc].index;

    if (cLocalMaxima[globalMaximumProc].value <= EPS) {
      status = SUCCESS;
      for (long k = 0; k < M; k++) {
        if (k == s) continue;
        bsp_put(k * N + t, &status, &status, 0, sizeof(STATUS));
      }
    } else {
      // Broadcast e to rest of the column
      for (long k = 0; k < M; k++) {
        if (k == s) continue;
        bsp_put(k * N + t, &e, &e, 0, sizeof(long));
      }
    }
  }
}

void findRowPhase1(long s, long t, long eModN, long locRows, flt *b) {
  if (t == 0 && eModN != t)
    bsp_put(s * N + eModN, b, b, 0, locRows * sizeof(flt));
}

void findRowPhase2(
        long s, long t, long eModN, long eDivN, long locRows, flt **A, flt *b, IndexValuePair *baLocalMinima) {
  if (eModN == t) { // Now find minimum
    long localIndex = -1;
    flt localMin = F_MAX;
    for (long i = 0; i < locRows; i++) {
      if (A[i][eDivN] > EPS) {
        flt ba = b[i] / A[i][eDivN];
        if (ba < localMin) {
          localIndex = i;
          localMin = ba;
        }
      }
    }

    baLocalMinima[s] = {.index=localIndex == -1 ? -1 : M * localIndex + s, .value=localMin};

    // Distribute to the rest of column if valid
    for (long k = 0; k < M; k++) {
      if (k == s) continue;
      baLocalMinima[k] = {.index=-1, .value=F_MAX};
      if (localIndex > -1)
        bsp_put(k * N + t, &baLocalMinima[s], baLocalMinima, s * sizeof(IndexValuePair), sizeof(IndexValuePair));
    }
  }
}

void findRowPhase3(
        long s, long t, long eModN, IndexValuePair *baLocalMinima, STATUS &status, long e, long &l) {
  if (eModN == t) {
    long globalMinimumProc = 0;
    for (long k = 1; k < M; k++) {
      if (baLocalMinima[k].value < baLocalMinima[globalMinimumProc].value)
        globalMinimumProc = k;
    }

    l = baLocalMinima[globalMinimumProc].index;
    if (l == -1) {
      status = UNBOUNDED;
      for (long k = 0; k < N; k++) {
        if (k == t) continue;
        bsp_put(s * N + k, &status, &status, 0, sizeof(STATUS));
      }
    } else {
      // Now share l with the rest of the row
      for (long k = 0; k < N; k++) {
        if (k == t) continue;
        bsp_put(s * N + k, &l, &l, 0, sizeof(long));
      }
    }
  }
}

void pivotPhase1(long s, long t, long eModN, long eDivN, long locRows, flt &ce, flt *c, flt *aie, flt **A) {
  if (eModN == t) {
    // Start two-phase broadcast of column e
    if (s == 0) {
      ce = c[eDivN];
      for (long k = 0; k < N; k++) {
        if (k == t) continue;
        bsp_put(0 * N + k, &ce, &ce, 0, sizeof(flt));
      }
    }
    for (long i = 0; i < locRows; i++)
      aie[i] = A[i][eDivN];
  }
  bsp_broadcast(aie, locRows, s * N + eModN, s * N + 0, s * N + t, 1, N, 0);
}

void pivotPhase2(long s, long t, long eModN, long locRows, flt *aie){
  bsp_broadcast(aie, locRows, s * N + eModN, s * N + 0, s * N + t, 1, N, 1);
}

void pivotPhase3(
        long s, long t, long lModM, long lDivM, long locCols, long eModN, long eDivN, flt *aie, flt *alj, flt **A,
        flt *bl, flt *b) {
  if (lModM == s) {
    long i = lDivM;
    for (long j = 0; j < locCols; j++) {
      if (eModN == t && j == eDivN)
        A[i][j] = 1 / aie[i];
      else
        A[i][j] /= aie[i];
    }
    if (t == 0)
      b[i] /= aie[i];

    // Distribute row to the rest of the column
    for (long j = 0; j < locCols; j++)
      alj[j] = A[i][j];
    if (t == 0) {
      *bl = b[i];
      for (long k = 0; k < M; k++) {
        if (k == s) continue;
        bsp_put(k * N + t, bl, bl, 0, sizeof(flt));
      }
    }
  }
  bsp_broadcast(alj, locCols, lModM * N + t, 0 * N + t, s * N + t, N, M, 0);
}

void pivotPhase4(long s, long t, flt *alj, long locCols, long lModM) {
  bsp_broadcast(alj, locCols, lModM * N + t, 0 * N + t, s * N + t, N, M, 1);
}

void pivotPhase5(
        long s, long t, long locRows, long locCols, long eDivN, long eModN, long lModM, long lDivM,
        flt **A, const flt *alj, const flt *aie, flt *b, flt bl, flt *c, flt ce, flt &v) {
  for (long i = 0; i < locRows; i++) {
    if (lModM == s && i == lDivM) continue;
    for (long j = 0; j < locCols; j++) {
      if (eModN == t && j == eDivN)
        A[i][j] = -A[i][j] * alj[j];
      else
        A[i][j] -= aie[i] * alj[j];
    }
  }

  if (t == 0) {
    for (long i = 0; i < locRows; i++) {
      if (lModM == s && i == lDivM) continue;
      b[i] -= aie[i] * bl;
    }
  }
  if (s == 0) {
    for (long j = 0; j < locCols; j++) {
      if (eModN == t && j == eDivN)
        c[j] = -ce * alj[j];
      else
        c[j] -= ce * alj[j];
    }
    if (t == 0)
      v += ce * bl;
  }
}

result simplex(long M, long N, long s, long t, long m, long n, flt **A, flt *c, flt *b,
               long *Basis, long *NonBasis) {
  long locRows = nloc(M, s, m);
  long locCols = nloc(N, t, n);
  flt v = 0;
  auto *cLocalMaxima = new IndexValuePair[N];
  auto *baLocalMinima = new IndexValuePair[M];
  long e = -1; // Entering index
  long l = -1; // Leaving index
  flt *aie = new flt[locRows];
  flt *alj = new flt[locCols];
  flt bl = 0;
  flt ce = 0;
  STATUS status = RUNNING;

  for (long j = 0; j < n; j++)
    NonBasis[j] = j;
  for (long i = 0; i < m; i++)
    Basis[i] = n + i;

  bsp_push_reg(aie, locRows * sizeof(flt));
  bsp_push_reg(alj, locCols * sizeof(flt));
  bsp_push_reg(b, locRows * sizeof(flt));
  bsp_push_reg(&e, sizeof(long));
  bsp_push_reg(&l, sizeof(long));
  bsp_push_reg(&bl, sizeof(flt));
  bsp_push_reg(&status, sizeof(STATUS));
  bsp_push_reg(&ce, sizeof(flt));
  bsp_push_reg(baLocalMinima, M * sizeof(IndexValuePair));
  bsp_push_reg(cLocalMaxima, N * sizeof(IndexValuePair));
  bsp_sync();
  if (PROFILING) stepFinished(0, s, t);

  long iterations = 0;

  while (true) {
    if (PRINT_TABLES)
      print_tableau(s, t, m, n, A, c, b, v, Basis, NonBasis);

    // Find some e maximizing c[e]
    e = 0; // global column index
    findColPhase1(s, t, locCols, c, cLocalMaxima);

    bsp_sync();
    if (PROFILING) stepFinished(1, s, t);

    findColPhase2(s, t, cLocalMaxima, e, status);

    bsp_sync();
    if (PROFILING) stepFinished(2, s, t);

    if (status == SUCCESS) break;

    long eModN = e % N;
    long eDivN = e / N;

    findRowPhase1(s, t, eModN, locRows, b);

    bsp_sync();
    if (PROFILING) stepFinished(3, s, t);

    findRowPhase2(s, t, eModN, eDivN, locRows, A, b, baLocalMinima);
    pivotPhase1(s, t, eModN, eDivN, locRows, ce, c, aie, A);

    bsp_sync();
    if (PROFILING) stepFinished(4, s, t);

    findRowPhase3(s, t, eModN, baLocalMinima, status, e, l);
    pivotPhase2(s, t, eModN, locRows, aie);

    bsp_sync();
    if (PROFILING) stepFinished(5, s, t);
    if (s == 0 && t == 0 && PRINT_TABLES) printf("x%li enters; x%li leaves.\n", NonBasis[e], Basis[l]);

    if (status == UNBOUNDED) {
      v = F_MAX;
      break;
    }
    swap(Basis[l], NonBasis[e]);

    // Now we compute row l
    long lModM = l % M;
    long lDivM = l / M;

    pivotPhase3(s, t, lModM, lDivM, locCols, eModN, eDivN, aie, alj, A, &bl, b);

    bsp_sync();
    if (PROFILING) stepFinished(6, s, t);

    pivotPhase4(s, t, alj, locCols, lModM);

    bsp_sync();
    if (PROFILING) stepFinished(7, s, t);

    pivotPhase5(s, t, locRows, locCols, eDivN, eModN, lModM, lDivM, A, alj, aie, b, bl, c, ce, v);
    iterations++;
    if (PROFILING) stepFinished(8, s, t);

    if (DEBUG) {
      for (long i = 0; i < locRows; i++) {
        if (b[i] < 0) {
          printf("Negative part of solution: b%li = %lf\n", Basis[s + M*i], b[i]);
        }
      }
      if (s == 0 && t == 0 && iterations % 100 == 0) {
        printf("k=%li, v=%lf\n", iterations, v);
      }
    }
  }

  bsp_pop_reg(aie);
  bsp_pop_reg(alj);
  bsp_pop_reg(b);
  bsp_pop_reg(baLocalMinima);
  bsp_pop_reg(&e);
  bsp_pop_reg(&l);
  bsp_pop_reg(&bl);
  bsp_pop_reg(&status);

  bsp_pop_reg(cLocalMaxima);
  bsp_pop_reg(&ce);
  delete[] cLocalMaxima;

  delete[] baLocalMinima;
  delete[] alj;
  delete[] aie;

  return {
          .z =  v,
          .iters = iterations,
          .status = status
  };
}

void easy_test_one_proc() {
  bsp_begin(1);
  long p = 1; /* p=M*N */
  long pid = bsp_pid();

  long s = pid % M;
  long t = pid / M;

  int n = 3;
  int m = 3;
  flt dA[] = {
          1, 1, 3,
          2, 2, 5,
          4, 1, 2
  };
  flt *A[] = {
          &dA[0],
          &dA[3],
          &dA[6]
  };
  flt b[] = {
          30,
          24,
          36
  };
  flt c[] = {3, 1, 2};
  flt v = 0;

  long *Basis = new long[m];
  long *NonBasis = new long[n];

  simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);

  bsp_end();
}

void easy_test_two_cols() {
  bsp_begin(2);
  long p = 2;
  long pid = bsp_pid();

  long s = pid % M;
  long t = pid / M;

  int n = 3;
  int m = 3;
  long *Basis = new long[m];
  long *NonBasis = new long[n];

  if (t == 0) {
    flt dA[] = {
            1, 3,
            2, 5,
            4, 2
    };
    flt *A[] = {
            &dA[0],
            &dA[2],
            &dA[4]
    };
    flt b[] = {
            30,
            24,
            36
    };
    flt c[] = {3, 2};
    flt v = 0;

    simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  } else {
    flt dA[] = {
            1,
            2,
            1
    };
    flt *A[] = {
            &dA[0],
            &dA[1],
            &dA[2]
    };
    flt *b;
    flt c[] = {1};
    flt v = 0;
    simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  }
  bsp_end();
}

void easy_test_two_rows() {
  bsp_begin(2);
  long p = 2;
  long pid = bsp_pid();

  long s = pid % M;
  long t = pid / M;

  int n = 3;
  int m = 3;

  long *Basis = new long[m];
  long *NonBasis = new long[n];

  if (s == 0) {
    flt dA[] = {
            1, 1, 3,
            4, 1, 2
    };
    flt *A[] = {
            &dA[0],
            &dA[3]
    };
    flt b[] = {
            30,
            36
    };
    flt c[] = {3, 1, 2};
    flt v = 0;
    simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  } else {
    flt dA[] = {
            2, 2, 5
    };
    flt *A[] = {
            &dA[0]
    };
    flt b[] = {
            24
    };
    flt *c;
    flt v = 0;
    simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  }
  bsp_end();
}

void inputMatrixSize(long &n, long &m) {
  bsp_push_reg(&m, sizeof(long));
  bsp_push_reg(&n, sizeof(long));
  bsp_push_reg(&M, sizeof(long));
  bsp_push_reg(&N, sizeof(long));
  bsp_sync();
  if (bsp_pid() == 0) {
    if (nArg <= 0) {
      printf("Please enter matrix size nxm:\n");
      if (scanf("%ldx%ld", &m, &n) != 2) bsp_abort("Entered number invalid!\n");
    } else n = m = nArg;
    for (long q = 0; q < bsp_nprocs(); q++) {
      bsp_put(q, &M, &M, 0, sizeof(long));
      bsp_put(q, &N, &N, 0, sizeof(long));
      bsp_put(q, &m, &m, 0, sizeof(long));
      bsp_put(q, &n, &n, 0, sizeof(long));
    }
    if (nArg <= 0) {
      printf("Linear Optimization of %ld by %ld matrix\n", m, n);
      printf("using the %ld by %ld cyclic distribution\n", M, N);
    } else printf(R"({"M": %li, "N": %li, "n": %li, )", M, N, n);
  }
  bsp_sync();
  bsp_pop_reg(&n);
  bsp_pop_reg(&m);
  bsp_pop_reg(&N);
  bsp_pop_reg(&M);
  bsp_sync();
}

void distribute_and_run(long n, long m, flt **gA, flt *gb, flt *gc) {
  long s = bsp_pid() / N;
  long t = bsp_pid() % N;

  long nlr = nloc(M, s, m); /* number of local rows */
  long nlc = nloc(N, t, n); /* number of local columns */
  flt **A = matallocd(nlr, nlc);
  flt *b = new flt[nlr];
  flt *c = new flt[nlc];
  bsp_push_reg(A[0], nlr * nlc * sizeof(flt));
  bsp_push_reg(b, nlr * nlc * sizeof(flt));
  bsp_push_reg(c, nlr * nlc * sizeof(flt));
  bsp_sync();

  if (s == 0 && t == 0) {
    if (nArg <= 0) printf("Now distributing the problem...\n");
    long maxLRows = ceil(((flt) m) / M);
    long maxLCols = ceil(((flt) n) / N);
    flt **lA = matallocd(M * N, maxLRows * maxLCols);
    flt **lb = matallocd(M, maxLRows);
    flt **lc = matallocd(N, maxLCols);
    for (long i = 0; i < m; i++) {
      for (long j = 0; j < n; j++) {
        lA[(i % M) * N + j % N][i / M * nloc(N, j % N, n) + j / N] = gA[i][j];
      }
      lb[i % M][i / M] = gb[i];
    }
    for (long j = 0; j < n; j++)
      lc[j % N][j / N] = gc[j];

    matfreed(gA);
    delete[] gb;
    delete[] gc;

    for (long i = 0; i < M; i++) {
      for (long j = 0; j < N; j++) {
        bsp_put(i * N + j, lA[i * N + j], A[0], 0, nloc(M, i, m) * nloc(N, j, n) * sizeof(flt));
      }
      bsp_put(i * N + 0, lb[i], b, 0, nloc(M, i, m) * sizeof(flt));
    }
    for (long j = 0; j < N; j++)
      bsp_put(0 * N + j, lc[j], c, 0, nloc(N, j, n) * sizeof(flt));

    bsp_sync();
    matfreed(lA);
    matfreed(lb);
    matfreed(lc);
  } else {
    bsp_sync();
  }

  bsp_pop_reg(A[0]);
  bsp_pop_reg(b);
  bsp_pop_reg(c);
  bsp_sync();
  long *Basis = new long[m];
  long *NonBasis = new long[n];

  if (s == 0 && t == 0 && nArg <= 0)
    printf("Start of Linear Optimization\n");
  bsp_sync();
  if (s == 0 && t == 0)
    lastTime = bsp_time();
  flt time0 = lastTime;
  result res = simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  bsp_sync();
  flt time1 = bsp_time();

  if (s == 0 && t == 0) {
    if (nArg <= 0) {
      printf("End of Linear Optimization\n");
      if (res.status == SUCCESS)
        printf("Success with optimal value %lf", res.z);
      else if (res.status == UNBOUNDED)
        printf("Unbounded (detected at value %lf)", res.z);
      printf(" using %li iterations.\n", res.iters);
      printf("This took only %.6lf seconds.\n", time1 - time0);
      if (PROFILING) {
        for (long i = 0; i < NUM_STEPS; i++) {
          printf("Step %li took %.6lf seconds.\n", i, times[i]);
        }
      }
    } else printf("\"t\": %lf, \"k\": %li, \"z\": %lf},\n", time1 - time0, res.iters, res.z);
  }

  matfreed(A);
  delete[] b;
  delete[] c;
  delete[] Basis;
  delete[] NonBasis;
}

void simplex_from_file() {
  bsp_begin(M * N);

  long m, n;

  bsp_push_reg(&n, sizeof(long));
  bsp_push_reg(&m, sizeof(long));
  bsp_push_reg(&M, sizeof(long));
  bsp_push_reg(&N, sizeof(long));
  bsp_sync();

  flt **A, *b, *c;
  if (bsp_pid() == 0) {
    if (Gargc < 2) bsp_abort("Invalid number of arguments!");
    std::string inputStr = std::string(Gargv[1]);
    std::ifstream fileshape = std::ifstream(inputStr + "-shape.csv");
    if (!fileshape) std::cout << "Could not open file: " + inputStr + "-shape.csv";
    std::string cell;
    std::getline(fileshape, cell);
    m = std::stol(cell);
    std::getline(fileshape, cell);
    n = std::stol(cell);
    fileshape.close();

    A = matallocd(m, n);
    std::ifstream fileA = std::ifstream(inputStr + "-A.csv");
    if (!fileA) std::cout << "Could not open file: " + inputStr + "-A.csv";
    for (long i = 0; i < m; i++) {
      std::string cell;
      for (long j = 0; j < n; j++) {
        if (!std::getline(fileA, cell, ',')) printf("Couldn't read A");
        A[i][j] = std::stod(cell);
      }
    }
    fileA.close();

    b = new flt[m];
    std::ifstream fileb = std::ifstream(inputStr + "-b.csv");
    if (!fileb) std::cout << "Could not open file: " + inputStr + "-b.csv";
    for (long i = 0; i < m; i++) {
      if (!std::getline(fileb, cell)) printf("Couldn't read b");
      b[i] = std::stod(cell);
    }
    fileb.close();

    c = new flt[n];
    std::ifstream filec = std::ifstream(inputStr + "-c.csv");
    if (!filec) std::cout << "Could not open file: " + inputStr + "-c.csv";
    for (long j = 0; j < n; j++) {
      std::getline(filec, cell);
      c[j] = std::stod(cell);
    }
    filec.close();
  }

  if (bsp_pid() == 0) {
    for (long k = 0; k < bsp_nprocs(); k++) {
      bsp_put(k, &m, &m, 0, sizeof(long));
      bsp_put(k, &n, &n, 0, sizeof(long));
      bsp_put(k, &M, &M, 0, sizeof(long));
      bsp_put(k, &N, &N, 0, sizeof(long));
    }
  }
  bsp_sync();

  bsp_pop_reg(&m);
  bsp_pop_reg(&n);
  bsp_pop_reg(&M);
  bsp_pop_reg(&N);

  distribute_and_run(n, m, A, b, c);

  bsp_end();
}

void simplex_from_rand() {
  bsp_begin(M * N);
  long n, m;
  inputMatrixSize(n, m);

  flt **gA, *gb, *gc;
  if (bsp_pid() == 0) {
    std::uniform_real_distribution <flt> unif(-1., 1.);
    std::default_random_engine re(12345);

    gA = matallocd(m, n);
    gb = new flt[m];
    gc = new flt[n];
    for (long i = 0; i < m; i++) {
      for (long j = 0; j < n; j++)
        gA[i][j] = unif(re);
      gb[i] = (unif(re) + 1) / 2;
    }
    for (long j = 0; j < n; j++)
      gc[j] = unif(re);
  }

  distribute_and_run(n, m, gA, gb, gc);

  bsp_end();
}

int main(int argc, char **argv) {
  Gargc = argc;
  Gargv = argv;
  bsp_init(simplex_from_rand, argc, argv);

  if (argc == 4) {
    M = std::stol(argv[1]);
    N = std::stol(argv[2]);
    nArg = std::stol(argv[3]);
  } else {
    printf("Please enter number of processor rows M:\n");
    if (scanf("%ld", &M) != 1) bsp_abort("Invalid input!");
    printf("Please enter number of processor columns N:\n");
    if (scanf("%ld", &N) != 1) bsp_abort("Invalid input!");
  }

  if (M * N > bsp_nprocs()) {
    if (argc != 4) printf("Sorry, only %u processors available.\n", bsp_nprocs());
    exit(EXIT_FAILURE);
  }

  simplex_from_rand();
  exit(EXIT_SUCCESS);
}
