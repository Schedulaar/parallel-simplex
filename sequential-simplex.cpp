#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>

int n, m;
int index (i, j) {
  return i*n + j;
}

void switchNumbers (int * x, int * y) {
  int tmp = *x;
  *x = *y;
  *y = tmp;
}

/**
 * e and l are indices of N and B.
 */
void pivot (
        int * N, int * B, double ** A, double * b, double * c, double v, int l, int e, // inputs
        int * nN, int * nB, double ** nA, double *nB, double * nc, double * nv         // outputs
      ) {
  double * nA = new double[n*m];
  double * nb = new double[m];
  double * nc = new double[m];
  double nv;
  int * nN = new int[m];
  int * nB = new int[n]

  // We will put x_N[e] into row l of the matrix and conversely put x_B[l] into column e.

  // Compute the coefficients of the equation for new basic variable x_N[e].
  nb[l] = b[l] / A[index(l,e)];
  for (int j = 0; j < m; j++) {
    if (j == e) continue;
    nA[index(l,j)] = A[index(l,j)] / A[index(l,e)];
  }
  nA[index(l,e)] = 1/A[index(l,e)];

  // Compute the coefficients of the remaining constraints.
  for (int i = 0; i < n; i++) {
    if (i == l) continue;
    nb[i] = b[i] - A[index(i,e)] * nb[l];
    for (int j = 0; j < m; j++) {
      if (j == e) continue;
      nA[index(i,j)] = A[index(i,j)] - A[index(i,e)]*nA[index(l,j)]
    }
    nA[index(i,e)] = -A[index(i,e)]*nA[index(l, e)];
  }

  // Compute the objective function
  *nv = v + c[N[e]] * nb[l];
  for (int j = 0; j < m; j++) {
    if (j == e) continue;
    nc[j] = c[j] - c[e]*nA[index(l,j)];
  }
  nc[e] = -c[e]*nA[index(l,e)];

  // Compute new sets of basic and nonbasic variables
  for (int i = 0; i < n; i++) {
    nB[i] = (i == l) ? N[e] : B[i];
  }
  for (int j = 0; j < m; j++) {
    nN[j] = (j == e) ? B[l] : N[j];
  }
}

int main(int argc, char **argv) {

  return 0;
}
