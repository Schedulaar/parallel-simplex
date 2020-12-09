#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>

bool PRINT = false;

int index(int i, int j, int rowLength) {
  return i * rowLength + j;
}

double EPS = 1e-20;

void swap(int &x, int &y) {
  int tmp = x;
  x = y;
  y = tmp;
}

void print_slack(int N[], int B[], double A[], double b[], double c[], double v, int n, int m) {
  printf("--------------------\n");
  printf("z  = %+.2f", v);
  for (int j = 0; j < n; j++)
    printf(" %+.2f*x%i", c[j], N[j]);
  printf("\n");
  for (int i = 0; i < m; i++) {
    printf("x%i = %+.2f", B[i], b[i]);
    for (int j = 0; j < n; j++)
      printf(" %+.2f*x%i", -A[index(i, j, n)], N[j]);
    printf("\n");
  }
}

/**
 * The variables e and l are indices of N and B.
 * It is possible to use the input variables as output variables as well.
 */
void pivot(
        int N[], int B[], double A[], double b[], double c[], double v, int l, int e, int n, int m, // inputs
        int oN[], int oB[], double oA[], double ob[], double oc[], double &ov                       // outputs
) {
  // We will put x_N[e] into row l of the matrix and conversely put x_B[l] into column e.

  // Compute the coefficients of the equation for new basic variable x_N[e].
  ob[l] = b[l] / A[index(l, e, n)];
  for (int j = 0; j < n; j++) {
    if (j == e) continue;
    oA[index(l, j, n)] = A[index(l, j, n)] / A[index(l, e, n)];
  }
  oA[index(l, e, n)] = 1 / A[index(l, e, n)];

  // Compute the coefficients of the remaining constraints.
  for (int i = 0; i < m; i++) {
    if (i == l) continue;
    ob[i] = b[i] - A[index(i, e, n)] * ob[l];
    for (int j = 0; j < n; j++) {
      if (j == e) continue;
      oA[index(i, j, n)] = A[index(i, j, n)] - A[index(i, e, n)] * oA[index(l, j, n)];
    }
    oA[index(i, e, n)] = -A[index(i, e, n)] * oA[index(l, e, n)];
  }

  // Compute the objective function
  ov = v + c[e] * ob[l];
  for (int j = 0; j < n; j++) {
    if (j == e) continue;
    oc[j] = c[j] - c[e] * oA[index(l, j, n)];
  }
  oc[e] = -c[e] * oA[index(l, e, n)];

  // Compute new sets of basic and nonbasic variables
  if (oB != B || N != oN) {
    for (int j = 0; j < m; j++) {
      oB[j] = B[j];
    }
    for (int i = 0; i < n; i++) {
      oN[i] = N[i];
    }
  }
  swap(oB[l], oN[e]);
}

std::string simplex_slack(
        int N[], int B[], double A[], double b[], double c[], double v, int n, int m, // inputs
        double x[], double &z                                                         // outputs
) {
  if (PRINT)
    print_slack(N, B, A, b, c, v, n, m);
  int e = 0;
  for (int j = 1; j < m; j++)
    e = c[j]>c[e] ? j : e;
  while (c[e] > EPS) { // While there exists j with c_j > 0

    // Find i in B minimizing b_i / a_ie
    int imin = -1;
    for (int i = 0; i < n; i++) {
      if (A[index(i, e, n)] > EPS &&
          (imin == -1 || b[imin] / A[index(imin, e, n)] > b[i] / A[index(i, e, n)]))
        imin = i;
    }
    if (imin == -1) return "unbounded";

    if (PRINT) printf("x%i enters; x%i leaves.\n", N[e], B[imin]);

    pivot(N, B, A, b, c, v, imin, e, n, m, // inputs
          N, B, A, b, c, v);            // outputs

    if (PRINT) print_slack(N, B, A, b, c, v, n, m);

    e = 0;
    for (int j = 1; j < m; j++)
      e = c[j]>c[e] ? j : e;
  }

  for (int j = 0; j < m; j++)
    x[B[j]] = b[j];
  for (int i = 0; i < n; i++)
    x[N[i]] = 0;
  z = v;

  return "success";
}

std::string simplex(
        double A[], double b[], double c[], int n, int m,
        double x[], double &z
) {
  int *N = new int[n];
  int *B = new int[m];
  double v;

  // PHASE 1: Find initial solution.

  // Find k minimizing b[i]
  int k = 0;
  for (int i = 1; i < m; i++)
    k = b[i] < b[k] ? i : k;
  if (b[k] >= 0) { // We can use the canonical slack form initially.
    for (int i = 0; i < n; i++)
      N[i] = i;
    for (int j = 0; j < m; j++)
      B[j] = n + j;
    v = 0;
  } else { // We have to solve an auxiliary problem, to find an initial solution.
    double *auxA = new double[m * (n + 1)];
    double *auxc = new double[n + 1];
    double *auxb = new double[m];
    int *auxN = new int[n + 1];
    int *auxB = new int[m];

    for (int i = 0; i < n + 1; i++) {
      auxc[i] = i == n ? -1. : 0.;
      auxN[i] = i == n ? n+m : i;
    }
    for (int j = 0; j < m; j++) {
      auxB[j] = n + j;
      auxb[j] = b[j];
    }
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n + 1; j++) {
        auxA[index(i, j, n + 1)] = j == n ? -1. : A[index(i, j, n)];
      }
    }

    double auxv = 0;
    pivot(auxN, auxB, auxA, auxb, auxc, auxv, k, n, n + 1, m,
          auxN, auxB, auxA, auxb, auxc, auxv);
    double auxz;
    double *xInit = new double[n + m + 1];
    std::string result = simplex_slack(auxN, auxB, auxA, auxb, auxc, auxv, n + 1, m,
                                       xInit, auxz);
    if (result != "success") {
      printf("Could not solve the auxiliary problem for finding an initial value.\n"
             "This should never happen.");
      return "error";
    }
    if (auxz > EPS) return "infeasible";

    // Check whether n+m is a basis variable. If so, do a non-degenerate pivot.
    for (int l = 0; l < m; l++) {
      if (auxB[l] == n + m){
        int e = 0;
        for (int j = 1; j < m; j++) {
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

    int remCol = n + 1;
    for (int j = 0; j < n + 1; j++) {
      if (auxN[j] == n+m)
        remCol = j;
      else
        N[j > remCol ? j - 1 : j] = auxN[j];
    }
    for (int i = 0; i < m; i++) {
      B[i] = auxB[i];
      b[i] = auxb[i];
    }

    // Remove column where variable x_{n+m} was in
    for (int j = 0; j < n + 1; j++) {
      if (j == remCol) continue;
      for (int i = 0; i < m; i++) {
        A[index(i, j - (j > remCol), n)] = auxA[index(i, j, n + 1)];
      }
    }

    // Build c
    double *newc = new double[n];
    for (int j = 0; j < n; j++) {
      newc[j] = N[j] < n ? c[N[j]] : 0;
      for (int i = 0; i < m; i++) {
        if (B[i] < n)
          newc[j] -= c[B[i]] * A[index(j, i, n)];
      }
    }

    // Build v
    v = 0;
    for (int i = 0; i < m; i++) {
      if (B[i] < n)
        v += c[B[i]] * b[i];
    }

    for (int j = 0; j < n; j++) c[j] = newc[j];
    delete[] newc;
  }

  // PHASE 2: Find optimal solution.
  return simplex_slack(N, B, A, b, c, v, n, m,
                       x, z);
}

void test() {
  int n = 3;
  int m = 3;
  // inputs
  int N[] = {0, 1, 2};
  int B[] = {3, 4, 5};
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
  int l = 2;
  int e = 0;

  // outputs
  double nA[n * m];
  double nb[m];
  double nc[m];
  double nv;
  int nN[m];
  int nB[n];

  double x[n + m];
  double z;
  // pivot(N, B, A, b, c, v, l, e, N, B, A, b, c, v);

  std::string result = simplex_slack(N, B, A, b, c, v, n, m,
                                     x, z);
  if (result == "success") {
    printf("Found optimal solution with objective value %lf:\n", z);
    for (int k = 0; k < n + m; k++) {
      printf("x%i=%lf,\n", k, x[k]);
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
  double *x = new double[4];
  double z;
  simplex(A, b, c, 2, 2,
          x, z);
  printf("%lf, %lf: %lf\n", x[0], x[1], z);
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

  double * x = new double[m+n];
  double z;

  printf("Starting Linear Optimization...\n");
  auto start = std::chrono::high_resolution_clock::now();

  simplex(A, b, c, n, m, x, z);

  auto finish = std::chrono::high_resolution_clock::now();
  printf("Finished Linear Optimization with optimal value %lf\n", z);
  std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count() << "ns\n";

  delete[] x;
  delete[] A;
  delete[] b;
  delete[] c;
}

int main(int argc, char **argv) {
  testFromFile();
  return 0;
}
