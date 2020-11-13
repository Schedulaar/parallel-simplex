#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>

int n, m;
int index (int i, int j) {
  return i*n + j;
}

void swap (int & x, int & y) {
  int tmp = x;
  x = y;
  y = tmp;
}

/**
 * e and l are indices of N and B.
 */
void pivot (
        int N[], int B[], double A[], double b[], double c[], double v, int l, int e, // inputs
        int nN[], int nB[], double nA[], double nb[], double nc[], double & nv         // outputs
      ) {
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
      nA[index(i,j)] = A[index(i,j)] - A[index(i,e)]*nA[index(l,j)];
    }
    nA[index(i,e)] = -A[index(i,e)]*nA[index(l, e)];
  }

  // Compute the objective function
  nv = v + c[e] * nb[l];
  for (int j = 0; j < m; j++) {
    if (j == e) continue;
    nc[j] = c[j] - c[e]*nA[index(l,j)];
  }
  nc[e] = -c[e]*nA[index(l,e)];

  // Compute new sets of basic and nonbasic variables
  for (int i = 0; i < n; i++) {
    nB[i] = B[i];
  }
  for (int j = 0; j < m; j++) {
    nN[j] = N[j];
  }
  swap(nB[l], nN[e]);
}

void test () {
  n = 3;
  m = 3;
  // inputs
  int N[] = { 0,1,2 };
  int B[] = { 3,4,5 };
  double A[] = {
          1,1,3,
          2,2,5,
          4,1,2
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
  double nA[n*m];
  double nb[m];
  double nc[m];
  double nv;
  int nN[m];
  int nB[n];

  pivot(N, B, A, b, c, v, l, e, N, B, A, b, c, v);
  printf("%i,%i,%i", N[0],c[1],v);
  printf("%i,%i,%i", A[0],B[1],b[2]);
}

int main(int argc, char **argv) {
  test();
  return 0;
}
