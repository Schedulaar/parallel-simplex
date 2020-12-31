#include <cstdio>
#include <cmath>
#include <string>
#include <bsp.hpp>
#include <random>
#include <iostream>
#include <fstream>

struct result {
  double z;
  long iters;
};

long M, N;
long nArg = -1;

enum STATUS {
  RUNNING, SUCCESS, UNBOUNDED
};
bool PRINT_TABLES = false;
bool PROFILING = false;
double lastTime;
double times[] = {0, 0, 0, 0, 0, 0, 0, 0};

void stepFinished(long step, long s, long t) {
  if (s == 0 && t == 0) {
    double newTime = bsp_time();
    times[step] += newTime - lastTime;
    lastTime = newTime;
  }
}

long min(long a, long b) { return a <= b ? a : b; }

void bsp_broadcast(double *x, long n, long src, long s0, long s, long stride, long p0, long phase) {
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
        bsp_put(dest, &x[t1 * b], x, t1 * b * sizeof(double), min(b, n - t1 * b) * sizeof(double));
    }
  }
}

double **matallocd(size_t m, size_t n) {
  /* This function allocates an m x n matrix of doubles */
  size_t i;
  double *pd, **ppd;
  ppd = new double *[m];
  if (ppd == NULL)
    bsp_abort("matallocd: not enough memory");
  pd = new double[m * n];
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

void matfreed(double **ppd) {
  /* This function frees a matrix of doubles */
  if (ppd != NULL) {
    if (ppd[0] != NULL)
      free(ppd[0]);
    free(ppd);
  }
}

void print_tableau(long M, long N, long s, long t, long m, long n, double **A, double *c, double *b, double v,
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
  }
}

result simplex(long M, long N, long s, long t, long m, long n, double **A, double *c, double *b,
               long *Basis, long *NonBasis) {
  long locRows = nloc(M, s, m);
  long locCols = nloc(N, t, n);
  double v = 0;
  double *cLocalMaxima;
  long *cLocalMaximaIndices;
  double *baLocalMinima = new double[M];
  long *baLocalMinimaIndices = new long[M];;
  long e = -1; // Entering index
  long l = -1; // Leaving index
  double *aie = new double[locRows];
  double *alj = new double[locCols];
  double bl;
  double ce;
  int status = RUNNING;
  double EPS = 1e-20;

  for (long j = 0; j < n; j++)
    NonBasis[j] = j;
  for (long i = 0; i < m; i++)
    Basis[i] = n + i;


  bsp_push_reg(aie, locRows * sizeof(double));
  bsp_push_reg(alj, locCols * sizeof(double));
  bsp_push_reg(b, locRows * sizeof(double));
  bsp_push_reg(baLocalMinima, M * sizeof(double));
  bsp_push_reg(baLocalMinimaIndices, M * sizeof(long));
  bsp_push_reg(&e, sizeof(long));
  bsp_push_reg(&l, sizeof(long));
  bsp_push_reg(&bl, sizeof(double));
  bsp_push_reg(&status, sizeof(int));
  cLocalMaxima = new double[N];
  cLocalMaximaIndices = new long[N];
  bsp_push_reg(cLocalMaxima, N * sizeof(double));
  bsp_push_reg(cLocalMaximaIndices, N * sizeof(long)); // global indices of local maxima
  bsp_push_reg(&ce, sizeof(double));
  bsp_sync();
  if (PROFILING) stepFinished(0, s, t);

  long iterations = 0;

  while (true) {
    if (PRINT_TABLES)
      print_tableau(M, N, s, t, m, n, A, c, b, v, Basis, NonBasis);

    // Find some e maximizing c[e]
    e = 0; // global column index
    if (s == 0) {
      cLocalMaxima[t] = c[0];
      cLocalMaximaIndices[t] = t;
      for (long j = 1; j < locCols; j++) {
        if (c[j] > cLocalMaxima[t]) {
          cLocalMaximaIndices[t] = t + N * j;
          cLocalMaxima[t] = c[j];
        }
      }
      // Broadcast local maximum
      for (long k = 0; k < N; k++) {
        if (k == t) continue;
        bsp_put(0 * N + k, &cLocalMaxima[t], cLocalMaxima, t * sizeof(double), sizeof(double));
        bsp_put(0 * N + k, &cLocalMaximaIndices[t], cLocalMaximaIndices, t * sizeof(long), sizeof(long));
      }
      bsp_sync();
      if (PROFILING) stepFinished(1, s, t);

      long globalMaximumProc = 0;
      for (long k = 1; k < N; k++) {
        if (cLocalMaxima[globalMaximumProc] < cLocalMaxima[k])
          globalMaximumProc = k;
      }
      e = cLocalMaximaIndices[globalMaximumProc];


      if (cLocalMaxima[globalMaximumProc] <= EPS) {
        status = SUCCESS;
        for (long k = 0; k < M; k++) {
          if (k == s) continue;
          bsp_put(k * N + t, &status, &status, 0, sizeof(STATUS));
        }
        if (PRINT_TABLES && t == 0)
          printf("Found an optimal solution!\n");
        bsp_sync();
        if (PROFILING) stepFinished(2, s, t);
        break;
      }

      // Broadcast e to rest of the column
      for (long k = 0; k < M; k++) {
        if (k == s) continue;
        bsp_put(k * N + t, &e, &e, 0, sizeof(long));
      }
      bsp_sync();
      if (PROFILING) stepFinished(2, s, t);
    } else {
      bsp_sync();
      if (PROFILING) stepFinished(1, s, t);
      bsp_sync();
      if (PROFILING) stepFinished(2, s, t);
      if (status == SUCCESS) break;
    }

    // Find l with a_el > 0 minimizing b_l / a_el
    // Do so by first distributing b_i to a_il
    // At the same time, we share column e!
    long eModN = e % N;
    long edivN = e / N;

    if (t == 0 && eModN != t)
      bsp_put(s * N + eModN, b, b, 0, locRows * sizeof(double));
    if (eModN == t) {
      // Start two-phase broadcast of column e
      if (s == 0) {
        ce = c[edivN];
        for (long k = 0; k < N; k++) {
          if (k == t) continue;
          bsp_put(0 * N + k, &ce, &ce, 0, sizeof(double));
        }
      }
      for (long i = 0; i < locRows; i++)
        aie[i] = A[i][edivN];
      bsp_sync();
      if (PROFILING) stepFinished(3, s, t);

      // Now find maximum
      baLocalMinima[s] = -1;
      baLocalMinimaIndices[s] = -1;
      for (long i = 0; i < locRows; i++) {
        if (A[i][edivN] > EPS) {
          double ba = b[i] / A[i][edivN];
          if (ba < baLocalMinima[s] || baLocalMinimaIndices[s] == -1) {
            baLocalMinima[s] = ba;
            baLocalMinimaIndices[s] = M * i + s;
          }
        }
      }

      // Distribute to the rest of column
      for (long k = 0; k < M; k++) {
        if (k == s) continue;
        bsp_put(k * N + t, &baLocalMinima[s], baLocalMinima, s * sizeof(double), sizeof(double));
        bsp_put(k * N + t, &baLocalMinimaIndices[s], baLocalMinimaIndices, s * sizeof(long), sizeof(long));
      }
      bsp_broadcast(aie, locRows, s * N + eModN, s * N + 0, s * N + t, 1, N, 0);
      bsp_sync();
      if (PROFILING) stepFinished(4, s, t);

      bsp_broadcast(aie, locRows, s * N + eModN, s * N + 0, s * N + t, 1, N, 1);

      long globalMinimumProc = -1;
      for (long k = 0; k < M; k++) {
        if (baLocalMinimaIndices[k] != -1 &&
            (globalMinimumProc == -1 || baLocalMinima[k] < baLocalMinima[globalMinimumProc]))
          globalMinimumProc = k;
      }
      if (globalMinimumProc == -1) {
        status = UNBOUNDED;
        for (long k = 0; k < N; k++) {
          if (k == t) continue;
          bsp_put(s * N + k, &status, &status, 0, sizeof(STATUS));
        }
        if (PRINT_TABLES && s == 0) {
          printf("Problem unbounded!\n");
        }
        bsp_sync();
        if (PROFILING) stepFinished(5, s, t);
        break;
      }
      l = baLocalMinimaIndices[globalMinimumProc];
      if (s == 0 && PRINT_TABLES) printf("%li enters; %li leaves.\n", e, l);

      // Now share l with the rest of the row
      for (long k = 0; k < N; k++) {
        if (k == t) continue;
        bsp_put(s * N + k, &l, &l, 0, sizeof(long));
      }
      bsp_sync();
      if (PROFILING) stepFinished(5, s, t);
    } else {
      bsp_sync();
      if (PROFILING) stepFinished(3, s, t);
      bsp_sync();
      if (PROFILING) stepFinished(4, s, t);

      bsp_broadcast(aie, locRows, s * N + eModN, s * N + 0, s * N + t, 1, N, 1);
      bsp_sync();
      if (PROFILING) stepFinished(5, s, t);

      if (status == UNBOUNDED) break;
    }
    swap(Basis[l], NonBasis[e]);

    // Now we compute row l
    long lModM = l % M;
    long lDivM = l / M;
    if (lModM == s) { // Then processor handles row l
      long i = lDivM;
      for (long j = 0; j < locCols; j++) {
        if (eModN == t && j == edivN)
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
        bl = b[i];
        for (long k = 0; k < M; k++) {
          if (k == s) continue;
          bsp_put(k * N + t, &bl, &bl, 0, sizeof(double));
        }
      }
    }

    bsp_broadcast(alj, locCols, lModM * N + t, 0 * N + t, s * N + t, N, M, 0);
    bsp_sync();
    if (PROFILING) stepFinished(6, s, t);

    bsp_broadcast(alj, locCols, lModM * N + t, 0 * N + t, s * N + t, N, M, 1);
    bsp_sync();
    if (PROFILING) stepFinished(7, s, t);

    // Compute rest of the constraints
    // Compute the coefficients of the remaining constraints.
    for (long i = 0; i < locRows; i++) {
      if (lModM == s && i == lDivM) continue;
      if (t == 0)
        b[i] -= aie[i] * bl;
      for (long j = 0; j < locCols; j++) {
        if (eModN == t && j == edivN)
          A[i][j] = -A[i][j] * alj[j];
        else
          A[i][j] -= aie[i] * alj[j];
      }
    }
    if (s == 0) {
      for (long j = 0; j < locCols; j++) {
        if (eModN == t && j == edivN)
          c[j] = -ce * alj[j];
        else
          c[j] -= ce * alj[j];
      }
      if (t == 0)
        v += ce * bl;
    }
    iterations++;
  }

  bsp_pop_reg(aie);
  bsp_pop_reg(alj);
  bsp_pop_reg(b);
  bsp_pop_reg(baLocalMinima);
  bsp_pop_reg(baLocalMinimaIndices);
  bsp_pop_reg(&e);
  bsp_pop_reg(&l);
  bsp_pop_reg(&bl);
  bsp_pop_reg(&status);

  if (s == 0) {
    bsp_pop_reg(cLocalMaxima);
    bsp_pop_reg(cLocalMaximaIndices);
    bsp_pop_reg(&ce);
    delete[] cLocalMaxima;
    delete[] cLocalMaximaIndices;
  }

  delete[] baLocalMinima;
  delete[] baLocalMinimaIndices;
  delete[] alj;
  delete[] aie;

  return {
          .z =  v,
          .iters = iterations
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
  double dA[] = {
          1, 1, 3,
          2, 2, 5,
          4, 1, 2
  };
  double *A[] = {
          &dA[0],
          &dA[3],
          &dA[6]
  };
  double b[] = {
          30,
          24,
          36
  };
  double c[] = {3, 1, 2};
  double v = 0;

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
    double dA[] = {
            1, 3,
            2, 5,
            4, 2
    };
    double *A[] = {
            &dA[0],
            &dA[2],
            &dA[4]
    };
    double b[] = {
            30,
            24,
            36
    };
    double c[] = {3, 2};
    double v = 0;

    simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  } else {
    double dA[] = {
            1,
            2,
            1
    };
    double *A[] = {
            &dA[0],
            &dA[1],
            &dA[2]
    };
    double *b;
    double c[] = {1};
    double v = 0;
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
    double dA[] = {
            1, 1, 3,
            4, 1, 2
    };
    double *A[] = {
            &dA[0],
            &dA[3]
    };
    double b[] = {
            30,
            36
    };
    double c[] = {3, 1, 2};
    double v = 0;
    simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  } else {
    double dA[] = {
            2, 2, 5
    };
    double *A[] = {
            &dA[0]
    };
    double b[] = {
            24
    };
    double *c;
    double v = 0;
    simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  }
  bsp_end();
}

long inputMatrixSize() {
  long n;
  bsp_push_reg(&n, sizeof(long));
  bsp_push_reg(&M, sizeof(long));
  bsp_push_reg(&N, sizeof(long));
  bsp_sync();
  if (bsp_pid() == 0) {
    if (nArg <= 0) {
      printf("Please enter matrix size n (nxn):\n");
      if (scanf("%ld", &n) != 1) bsp_abort("Entered number invalid!\n");
    } else n = nArg;
    for (long q = 0; q < bsp_nprocs(); q++) {
      bsp_put(q, &M, &M, 0, sizeof(long));
      bsp_put(q, &N, &N, 0, sizeof(long));
      bsp_put(q, &n, &n, 0, sizeof(long));
    }
    if (nArg <= 0) {
      printf("Linear Optimization of %ld by %ld matrix\n", n, n);
      printf("using the %ld by %ld cyclic distribution\n", M, N);
    } else printf(R"({"M": %li, "N": %li, "n": %li, )", M, N, n);
  }
  bsp_sync();
  bsp_pop_reg(&n);
  bsp_pop_reg(&N);
  bsp_pop_reg(&M);
  bsp_sync();
  return n;
}

void distribute_and_run(long n, long m, double **gA, double *gb, double *gc) {
  long s = bsp_pid() / N;
  long t = bsp_pid() % N;

  long nlr = nloc(M, s, m); /* number of local rows */
  long nlc = nloc(N, t, n); /* number of local columns */
  double **A = matallocd(nlr, nlc);
  double *b = new double[nlr];
  double *c = new double[nlc];
  bsp_push_reg(A[0], nlr * nlc * sizeof(double));
  bsp_push_reg(b, nlr * nlc * sizeof(double));
  bsp_push_reg(c, nlr * nlc * sizeof(double));
  bsp_sync();

  if (s == 0 && t == 0) {
    if (nArg <= 0) printf("Now distributing the problem...\n");
    long maxLRows = ceil(((double) m) / M);
    long maxLCols = ceil(((double) n) / N);
    double **lA = matallocd(M * N, maxLRows * maxLCols);
    double **lb = matallocd(M, maxLRows);
    double **lc = matallocd(N, maxLCols);
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
        bsp_put(i * N + j, lA[i * N + j], A[0], 0, nloc(M, i, m) * nloc(N, j, n) * sizeof(double));
      }
      bsp_put(i * N + 0, lb[i], b, 0, nloc(M, i, m) * sizeof(double));
    }
    for (long j = 0; j < N; j++)
      bsp_put(0 * N + j, lc[j], c, 0, nloc(N, j, n) * sizeof(double));

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
  double time0 = lastTime;
  result res = simplex(M, N, s, t, m, n, A, c, b, Basis, NonBasis);
  bsp_sync();
  double time1 = bsp_time();

  if (s == 0 && t == 0) {
    if (nArg <= 0) {
      printf("End of Linear Optimization: Optimal value %lf using %li iterations.\n", res.z, res.iters);
      printf("This took only %.6lf seconds.\n", time1 - time0);
      if (PROFILING) {
        for (long i = 0; i < 8; i++) {
          printf("Step %li took %.6lf seconds.\n", i, times[i]);
        }
      }
    } else printf("\"t\": %lf},\n", time1 - time0);
  }

  matfreed(A);
  delete[] b;
  delete[] c;
  delete[] Basis;
  delete[] NonBasis;
}

void simplex_from_file() {
  bsp_begin(M * N);

  long n, m;
  n = m = inputMatrixSize();

  double **gA, *gb, *gc;
  if (bsp_pid() == 0) {
    long maxLRows = ceil(((double) m) / M);
    long maxLCols = ceil(((double) n) / N);

    gA = matallocd(M * N, maxLRows * maxLCols);
    gA = matallocd(m, n);
    std::ifstream fileA;
    fileA.open(std::to_string(n) + "-A.csv");
    if (!fileA) bsp_abort(("Could not open file: " + std::to_string(n) + "-A.csv").c_str());
    for (long i = 0; i < m; i++) {
      std::string cell;
      for (long j = 0; j < n; j++) {
        if (!std::getline(fileA, cell, ',')) bsp_abort("Couldn't read A");
        gA[i][j] = std::stod(cell);
      }
    }
    fileA.close();


    gb = new double[m];
    std::ifstream fileb;
    fileb.open(std::to_string(n) + "-b.csv");
    if (!fileb) bsp_abort(("Could not open file: " + std::to_string(n) + "-b.csv").c_str());
    for (long i = 0; i < m; i++) {
      std::string cell;
      if (!std::getline(fileb, cell)) bsp_abort("Couldn't read b");
      gb[i] = std::stod(cell);
    }
    fileb.close();


    gc = new double[n];
    std::ifstream filec(std::to_string(n) + "-c.csv");
    for (long j = 0; j < n; j++) {
      std::string cell;
      if (!std::getline(filec, cell)) bsp_abort("Couldn't read c");
      gc[j] = std::stod(cell);
    }
    filec.close();
  }

  distribute_and_run(n, m, gA, gb, gc);

  bsp_end();
}

void simplex_from_rand() {
  bsp_begin(M * N);
  long n, m;
  n = m = inputMatrixSize();

  double **gA, *gb, *gc;
  if (bsp_pid() == 0) {
    srand(1);
    std::uniform_real_distribution<double> unif(0., 1.);
    std::default_random_engine re;

    gA = matallocd(m, n);
    gb = new double[m];
    gc = new double[n];
    for (long i = 0; i < m; i++) {
      for (long j = 0; j < n; j++)
        gA[i][j] = unif(re);
      gb[i] = unif(re);
    }
    for (long j = 0; j < n; j++)
      gc[j] = unif(re);
  }

  distribute_and_run(n, m, gA, gb, gc);

  bsp_end();
}

int main(int argc, char **argv) {
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
