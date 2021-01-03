#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>
#include <random>
#include <fstream>
#include <chrono>

using flt = float;

bool PRINT = false;

long index(long i, long j, long rowLength) {
  return i * rowLength + j;
}

flt EPS = 1e-20;

void swap(long &x, long &y) {
  long tmp = x;
  x = y;
  y = tmp;
}

void print_slack(long N[], long B[], flt A[], flt b[], flt c[], flt v, long n, long m) {
  printf("--------------------\n");
  printf("z  = %+.2f", v);
  for (long j = 0; j < n; j++)
    printf(" %+.2f*x%li", c[j], N[j]);
  printf("\n");
  for (long i = 0; i < m; i++) {
    printf("x%li = %+.2f", B[i], b[i]);
    for (long j = 0; j < n; j++)
      printf(" %+.2f*x%li", -A[index(i, j, n)], N[j]);
    printf("\n");
  }
}

/**
 * The variables e and l are indices of N and B.
 * It is possible to use the input variables as output variables as well.
 */
void pivot(
        long N[], long B[], flt A[], flt b[], flt c[], flt v, long l, long e, long n, long m, // inputs
        long oN[], long oB[], flt oA[], flt ob[], flt oc[], flt &ov                       // outputs
) {
  // We will put x_N[e] into row l of the matrix and conversely put x_B[l] into column e.

  // Compute the coefficients of the equation for new basic variable x_N[e].
  ob[l] = b[l] / A[index(l, e, n)];
  for (long j = 0; j < n; j++) {
    if (j == e) continue;
    oA[index(l, j, n)] = A[index(l, j, n)] / A[index(l, e, n)];
  }
  oA[index(l, e, n)] = 1 / A[index(l, e, n)];

  // Compute the coefficients of the remaining constraints.
  for (long i = 0; i < m; i++) {
    if (i == l) continue;
    ob[i] = b[i] - A[index(i, e, n)] * ob[l];
    for (long j = 0; j < n; j++) {
      if (j == e) continue;
      oA[index(i, j, n)] = A[index(i, j, n)] - A[index(i, e, n)] * oA[index(l, j, n)];
    }
    oA[index(i, e, n)] = -A[index(i, e, n)] * oA[index(l, e, n)];
  }

  // Compute the objective function
  ov = v + c[e] * ob[l];
  for (long j = 0; j < n; j++) {
    if (j == e) continue;
    oc[j] = c[j] - c[e] * oA[index(l, j, n)];
  }
  oc[e] = -c[e] * oA[index(l, e, n)];

  // Compute new sets of basic and nonbasic variables
  if (oB != B || N != oN) {
    for (long j = 0; j < m; j++) {
      oB[j] = B[j];
    }
    for (long i = 0; i < n; i++) {
      oN[i] = N[i];
    }
  }
  swap(oB[l], oN[e]);
}

std::string simplex_slack(
        long N[], long B[], flt A[], flt b[], flt c[], flt v, long n, long m, // inputs
        flt &z, long &iters                                                            // outputs
) {
  if (PRINT)
    print_slack(N, B, A, b, c, v, n, m);
  long e = 0;
  for (long j = 1; j < m; j++)
    e = c[j]>c[e] ? j : e;
  while (c[e] > EPS) { // While there exists j with c_j > 0

    // Find i in B minimizing b_i / a_ie
    long imin = -1;
    for (long i = 0; i < n; i++) {
      if (A[index(i, e, n)] > EPS &&
          (imin == -1 || b[imin] / A[index(imin, e, n)] > b[i] / A[index(i, e, n)]))
        imin = i;
    }
    if (imin == -1) return "unbounded";

    if (PRINT) printf("x%li enters; x%li leaves.\n", N[e], B[imin]);

    pivot(N, B, A, b, c, v, imin, e, n, m, // inputs
          N, B, A, b, c, v);            // outputs

    if (PRINT) print_slack(N, B, A, b, c, v, n, m);

    e = 0;
    for (long j = 1; j < m; j++)
      e = c[j]>c[e] ? j : e;
    iters++;
  }
  z = v;

  return "success";
}

std::string simplex(
        flt A[], flt b[], flt c[], long n, long m,
        flt &z, long* B, long &iters
) {
  long *N = new long[n];
  flt v;
  iters = 0;

  // PHASE 1: Find initial solution.

  // Find k minimizing b[i]
  long k = 0;
  for (long i = 1; i < m; i++)
    k = b[i] < b[k] ? i : k;
  if (b[k] >= 0) { // We can use the canonical slack form initially.
    for (long i = 0; i < n; i++)
      N[i] = i;
    for (long j = 0; j < m; j++)
      B[j] = n + j;
    v = 0;
  } else { // We have to solve an auxiliary problem, to find an initial solution.
    flt *auxA = new flt[m * (n + 1)];
    flt *auxc = new flt[n + 1];
    flt *auxb = new flt[m];
    long *auxN = new long[n + 1];
    long *auxB = new long[m];

    for (long i = 0; i < n + 1; i++) {
      auxc[i] = i == n ? -1. : 0.;
      auxN[i] = i == n ? n+m : i;
    }
    for (long j = 0; j < m; j++) {
      auxB[j] = n + j;
      auxb[j] = b[j];
    }
    for (long i = 0; i < m; i++) {
      for (long j = 0; j < n + 1; j++) {
        auxA[index(i, j, n + 1)] = j == n ? -1. : A[index(i, j, n)];
      }
    }

    flt auxv = 0;
    pivot(auxN, auxB, auxA, auxb, auxc, auxv, k, n, n + 1, m,
          auxN, auxB, auxA, auxb, auxc, auxv);
    flt auxz;
    std::string result = simplex_slack(auxN, auxB, auxA, auxb, auxc, auxv, n + 1, m,
                                       auxz, iters);
    if (result != "success") {
      printf("Could not solve the auxiliary problem for finding an initial value.\n"
             "This should never happen.");
      return "error";
    }
    if (auxz > EPS) return "infeasible";

    // Check whether n+m is a basis variable. If so, do a non-degenerate pivot.
    for (long l = 0; l < m; l++) {
      if (auxB[l] == n + m){
        long e = 0;
        for (long j = 1; j < m; j++) {
          if (std::abs(auxA[index(e, j, n + 1)]) > std::abs(auxA[index(e, l, n + 1)]))
            e = j;
        }
        pivot(auxN, auxB, auxA, auxb, auxc, auxv, l, e, n + 1, m,
              auxN, auxB, auxA, auxb, auxc, auxv);
        iters++;
        break;
      }
    }

    // Finally we have an initial solution.
    // Now, we remove variable n+m again.

    long remCol = n + 1;
    for (long j = 0; j < n + 1; j++) {
      if (auxN[j] == n+m)
        remCol = j;
      else
        N[j > remCol ? j - 1 : j] = auxN[j];
    }
    for (long i = 0; i < m; i++) {
      B[i] = auxB[i];
      b[i] = auxb[i];
    }

    // Remove column where variable x_{n+m} was in
    for (long j = 0; j < n + 1; j++) {
      if (j == remCol) continue;
      for (long i = 0; i < m; i++) {
        A[index(i, j - (j > remCol), n)] = auxA[index(i, j, n + 1)];
      }
    }

    // Build c
    flt *newc = new flt[n];
    for (long j = 0; j < n; j++) {
      newc[j] = N[j] < n ? c[N[j]] : 0;
      for (long i = 0; i < m; i++) {
        if (B[i] < n)
          newc[j] -= c[B[i]] * A[index(j, i, n)];
      }
    }

    // Build v
    v = 0;
    for (long i = 0; i < m; i++) {
      if (B[i] < n)
        v += c[B[i]] * b[i];
    }

    for (long j = 0; j < n; j++) c[j] = newc[j];
    delete[] newc;
  }

  // PHASE 2: Find optimal solution.
  return simplex_slack(N, B, A, b, c, v, n, m,
                       z, iters);
}

void test() {
  long n = 3;
  long m = 3;
  // inputs
  long N[] = {0, 1, 2};
  long B[] = {3, 4, 5};
  flt A[] = {
          1, 1, 3,
          2, 2, 5,
          4, 1, 2
  };
  flt b[] = {
          30,
          24,
          36
  };
  flt c[] = {3, 1, 2};
  flt v = 0;
  long l = 2;
  long e = 0;

  // outputs
  flt nA[n * m];
  flt nb[m];
  flt nc[m];
  flt nv;
  long nN[m];
  long nB[n];

  flt z;
  long iters = 0;
  // pivot(N, B, A, b, c, v, l, e, N, B, A, b, c, v);

  std::string result = simplex_slack(N, B, A, b, c, v, n, m,
                                     z, iters);
  if (result == "success") {
    printf("Found optimal solution with objective value %lf:\n", z);
    for (long k = 0; k < m; k++) {
      printf("x%li=%lf,\n", B[k], b[k]);
    }
  } else if (result == "unbounded") {
    printf("The problem is unbounded!\n");
  } else
    printf("Something definitely went wrong!\n");
}

void testPhase1() {
  flt A[] = {
          2, -1,
          1, -5
  };
  flt b[] = {
          2,
          -4
  };
  flt c[] = {2, -1};
  flt z;
  long B[2];
  long iters;
  simplex(A, b, c, 2, 2,
          z, B, iters);
  printf("Found optimal solution with objective value %lf:\n", z);
  for (long k = 0; k < 2; k++) {
    printf("x%li=%lf,\n", B[k], b[k]);
  }
}

void testFromFile() {
  long n;
  printf("Choose matrix size nxn!\n");
  if (scanf("%li", &n) != 1) printf("Invalid input!");
  long m = n;

  flt * A = new flt[m*n];
  std::ifstream fileA;
  fileA.open(std::to_string(n) + "-A.csv");
  if (!fileA) std::cout << "Could not open file: " + std::to_string(n) + "-A.csv";
  for (long i = 0; i < m; i++) {
    std::string cell;
    for (long j = 0; j < n; j++) {
      if(!std::getline(fileA, cell, ',')) printf("Couldn't read A");
      A[index(i,j,n)] = std::stod(cell);
    }
  }
  fileA.close();

  flt * b = new flt[m];
  std::ifstream fileb;
  fileb.open(std::to_string(n) + "-b.csv");
  if (!fileb) std::cout << "Could not open file: " + std::to_string(n) + "-b.csv";
  std::string cell;
  for (long i = 0; i < m; i++) {
    if(!std::getline(fileb, cell)) printf("Couldn't read b");
    b[i] = std::stod(cell);
  }
  fileb.close();

  flt * c = new flt[n];
  std::ifstream filec (std::to_string(n) + "-c.csv");
  for (long j = 0; j < n; j++) {
    std::string cell;
    std::getline(filec, cell);
    c[j] = std::stod(cell);
  }
  filec.close();

  flt z;
  long * B = new long[m];

  printf("Starting Linear Optimization...\n");
  long iters = 0;
  auto start = std::chrono::high_resolution_clock::now();

  simplex(A, b, c, n, m, z, B, iters);

  auto finish = std::chrono::high_resolution_clock::now();
  printf("Finished Linear Optimization with optimal value %lf\n", z);
  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << "ns\n";

  delete[] A;
  delete[] B;
  delete[] b;
  delete[] c;
}

void testFromRand (int argc, char ** argv) {
  long n;
  if (argc != 2) {
    printf("Choose matrix size nxn!\n");
    if (scanf("%li", &n) != 1) printf("Invalid input!");
  } else n = std::stol(argv[1]);
  long m = n;

  flt z;
  long * B = new long[m];

  std::uniform_real_distribution<flt> unif(0., 1.);
  std::default_random_engine re(12345);

  flt * A = new flt[m*n];
  flt * b = new flt[m];
  flt * c = new flt[n];
  for (long i = 0; i < m; i++) {
    for (long j = 0; j < n; j++)
      A[index(i,j,n)] = unif(re);
    b[i] = unif(re);
  }
  for (long j = 0; j < n; j++)
    c[j] = unif(re);


  if (argc != 2) printf("Starting Linear Optimization...\n");
  long iters = 0;
  auto start = std::chrono::high_resolution_clock::now();

  simplex(A, b, c, n, m, z, B, iters);

  auto finish = std::chrono::high_resolution_clock::now();
  if (argc != 2) printf("Finished Linear Optimization with optimal value %lf\n", z);
  double ns = std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();
  double s = ((double) ns) * 1e-9;
  if (argc != 2) printf("This took %lf seconds and %li iterations\n", s, iters);
  else printf(R"({"n": %li, "t": %lf, "k": %li, "z": %lf},)" "\n", n, s, iters, z);

  delete[] A;
  delete[] B;
  delete[] b;
  delete[] c;
}

int main(int argc, char **argv) {
  testFromRand(argc, argv);
  return 0;
}
