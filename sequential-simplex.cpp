#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>
#include <random>
#include <fstream>
#include <chrono>

bool PRINT = false;

long index(long i, long j, long rowLength) {
  return i * rowLength + j;
}

double EPS = 1e-20;

void swap(long &x, long &y) {
  long tmp = x;
  x = y;
  y = tmp;
}

void print_slack(long N[], long B[], double A[], double b[], double c[], double v, long n, long m) {
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
        long N[], long B[], double A[], double b[], double c[], double v, long l, long e, long n, long m, // inputs
        long oN[], long oB[], double oA[], double ob[], double oc[], double &ov                       // outputs
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
        long N[], long B[], double A[], double b[], double c[], double v, long n, long m, // inputs
        double &z                                                         // outputs
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
  }
  z = v;

  return "success";
}

std::string simplex(
        double A[], double b[], double c[], long n, long m,
        double &z, long* B
) {
  long *N = new long[n];
  double v;

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
    double *auxA = new double[m * (n + 1)];
    double *auxc = new double[n + 1];
    double *auxb = new double[m];
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

    double auxv = 0;
    pivot(auxN, auxB, auxA, auxb, auxc, auxv, k, n, n + 1, m,
          auxN, auxB, auxA, auxb, auxc, auxv);
    double auxz;
    std::string result = simplex_slack(auxN, auxB, auxA, auxb, auxc, auxv, n + 1, m,
                                       auxz);
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
    double *newc = new double[n];
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
                       z);
}

void test() {
  long n = 3;
  long m = 3;
  // inputs
  long N[] = {0, 1, 2};
  long B[] = {3, 4, 5};
  double A[] = {
          1, 1, 3,
          2, 2, 5,
          4, 1, 2
  };
  double b[] = {
          30,
          24,
          36
  };
  double c[] = {3, 1, 2};
  double v = 0;
  long l = 2;
  long e = 0;

  // outputs
  double nA[n * m];
  double nb[m];
  double nc[m];
  double nv;
  long nN[m];
  long nB[n];

  double z;
  // pivot(N, B, A, b, c, v, l, e, N, B, A, b, c, v);

  std::string result = simplex_slack(N, B, A, b, c, v, n, m,
                                     z);
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
  double A[] = {
          2, -1,
          1, -5
  };
  double b[] = {
          2,
          -4
  };
  double c[] = {2, -1};
  double z;
  long B[2];
  simplex(A, b, c, 2, 2,
          z, B);
  printf("Found optimal solution with objective value %lf:\n", z);
  for (long k = 0; k < 2; k++) {
    printf("x%li=%lf,\n", B[k], b[k]);
  }
}

void testFromFile() {
  long n;
  printf("Choose matrix size nxn!\n");
  scanf("%li", &n);
  long m = n;

  double * A = new double[m*n];
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

  double * b = new double[m];
  std::ifstream fileb;
  fileb.open(std::to_string(n) + "-b.csv");
  if (!fileb) std::cout << "Could not open file: " + std::to_string(n) + "-b.csv";
  std::string cell;
  for (long i = 0; i < m; i++) {
    if(!std::getline(fileb, cell)) printf("Couldn't read b");
    b[i] = std::stod(cell);
  }
  fileb.close();

  double * c = new double[n];
  std::ifstream filec (std::to_string(n) + "-c.csv");
  for (long j = 0; j < n; j++) {
    std::string cell;
    std::getline(filec, cell);
    c[j] = std::stod(cell);
  }
  filec.close();

  double z;
  long * B = new long[m];

  printf("Starting Linear Optimization...\n");
  auto start = std::chrono::high_resolution_clock::now();

  simplex(A, b, c, n, m, z, B);

  auto finish = std::chrono::high_resolution_clock::now();
  printf("Finished Linear Optimization with optimal value %lf\n", z);
  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << "ns\n";

  delete[] A;
  delete[] B;
  delete[] b;
  delete[] c;
}

void testFromRand () {
  long n;
  printf("Choose matrix size nxn!\n");
  scanf("%li", &n);
  long m = n;

  double z;
  long * B = new long[m];

  srand(1);
  std::uniform_real_distribution<double> unif(0., 1.);
  std::default_random_engine re;

  double * A = new double[m*n];
  double * b = new double[m];
  double * c = new double[n];
  for (long i = 0; i < m; i++) {
    for (long j = 0; j < n; j++)
      A[index(i,j,n)] = unif(re);
    b[i] = unif(re);
  }
  for (long j = 0; j < n; j++)
    c[j] = unif(re);


  printf("Starting Linear Optimization...\n");
  auto start = std::chrono::high_resolution_clock::now();

  simplex(A, b, c, n, m, z, B);

  auto finish = std::chrono::high_resolution_clock::now();
  printf("Finished Linear Optimization with optimal value %lf\n", z);
  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << "ns\n";

  delete[] A;
  delete[] B;
  delete[] b;
  delete[] c;
}

int main(int argc, char **argv) {
  testFromRand();
  return 0;
}
