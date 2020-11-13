#include <cstdio>
#include <string>

int index(int i, int j, int rowLength) {
  return i * rowLength + j;
}

double EPS = 1e-10;

void swap(int &x, int &y) {
  int tmp = x;
  x = y;
  y = tmp;
}

/**
 * The variables e and l are indices of N and B.
 * It is possible to use the input variables as output variables as well.
 */
void pivot(
        int N[], int B[], double A[], double b[], double c[], double v, int l, int e, int n, int m, // inputs
        int nN[], int nB[], double nA[], double nb[], double nc[], double &nv                       // outputs
) {
  // We will put x_N[e] into row l of the matrix and conversely put x_B[l] into column e.

  // Compute the coefficients of the equation for new basic variable x_N[e].
  nb[l] = b[l] / A[index(l, e, n)];
  for (int j = 0; j < n; j++) {
    if (j == e) continue;
    nA[index(l, j, n)] = A[index(l, j, n)] / A[index(l, e, n)];
  }
  nA[index(l, e, n)] = 1 / A[index(l, e, n)];

  // Compute the coefficients of the remaining constraints.
  for (int i = 0; i < m; i++) {
    if (i == l) continue;
    nb[i] = b[i] - A[index(i, e, n)] * nb[l];
    for (int j = 0; j < n; j++) {
      if (j == e) continue;
      nA[index(i, j, n)] = A[index(i, j, n)] - A[index(i, e, n)] * nA[index(l, j, n)];
    }
    nA[index(i, e, n)] = -A[index(i, e, n)] * nA[index(l, e, n)];
  }

  // Compute the objective function
  nv = v + c[e] * nb[l];
  for (int j = 0; j < n; j++) {
    if (j == e) continue;
    nc[j] = c[j] - c[e] * nA[index(l, j, n)];
  }
  nc[e] = -c[e] * nA[index(l, e, n)];

  // Compute new sets of basic and nonbasic variables
  if (nB != B || N != nN) {
    for (int j = 0; j < m; j++) {
      nB[j] = B[j];
    }
    for (int i = 0; i < n; i++) {
      nN[i] = N[i];
    }
  }
  swap(nB[l], nN[e]);
}

std::string simplex_slack(
        int N[], int B[], double A[], double b[], double c[], double v, int n, int m, // inputs
        double x[], double &z                                                         // outputs
) {
  int e;
  for (e = 0; e < m; e++)
    if (c[e] > EPS) break;
  while (e < m) { // While there exists j with c_j > 0

    // Find i in B minimizing b_i / a_ie
    int imin = -1;
    for (int i = 0; i < n; i++) {
      if (A[index(i, e, n)] > EPS &&
          (imin == -1 || b[imin] / A[index(imin, e, n)] > b[i] / A[index(i, e, n)]))
        imin = i;
    }
    if (imin == -1) return "unbounded";

    pivot(N, B, A, b, c, v, imin, e, n, m, // inputs
          N, B, A, b, c, v);            // outputs

    for (e = 0; e < m; e++)
      if (c[e] > EPS) break;
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
      auxN[i] = i;
    }
    for (int j = 0; j < m; j++) {
      auxB[j] = n + 1 + j;
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
    double *xInit = new double[n + 1 + m];
    std::string result = simplex_slack(auxN, auxB, auxA, auxb, auxc, auxv, n + 1, m,
                                       xInit, auxz);
    if (result != "success") {
      printf("Could not solve the auxiliary problem for finding an initial value.\n"
             "This should never happen.");
      return "error";
    }
    if (auxz > EPS) return "infeasible";

    // Check whether n is a basis variable. If so, do a non-degenerate pivot.
    int l = -1;
    for (int i = 0; i < m; i++) {
      if (auxB[i] != n) continue;
      l = i;
      break;
    }
    if (l != -1) {
      int e = 0;
      for (int j = 1; j < m; j++) {
        if (std::abs(auxA[index(e, j, n + 1)]) > std::abs(auxA[index(e, l, n + 1)]))
          e = j;
      }
      pivot(auxN, auxB, auxA, auxb, auxc, auxv, l, e, n + 1, m,
            auxN, auxB, auxA, auxb, auxc, auxv);
    }

    // Finally we have an initial solution.
    // Now, we remove variable n again.

    int remCol = n + 1;
    for (int j = 0; j < n + 1; j++) {
      if (auxN[j] == n)
        remCol = j;
      else
        N[j > remCol ? j - 1 : j] = auxN[j] < n ? auxN[j] : auxN[j] - 1;
    }
    for (int i = 0; i < m; i++) {
      B[i] = auxB[i] < n ? auxB[i] : auxB[i] - 1;
      b[i] = auxb[i];
    }

    // Remove column where variable x_n was in
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
    printf("Found optimal solution with object value %lf:\n", z);
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

int main(int argc, char **argv) {
  testPhase1();
  return 0;
}
