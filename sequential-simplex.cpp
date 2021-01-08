#include <cstdio>
#include <cmath>
#include <string>
#include <iostream>
#include <random>
#include <fstream>
#include <chrono>
#include <limits>
#include <string.h>

using flt = double;
flt INF = std::numeric_limits<flt>::infinity();
flt EPS = 1e-9;

bool PRINT = false;
bool DEBUG = false;

long index(long i, long j, long rowLength) {
  return i * rowLength + j;
}

flt max(flt x, flt y) { return x >= y ? x : y; }

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
        long N[], long B[], flt A[], flt b[], flt c[], flt v, flt c2[], flt v2, long l, long e, long n, long m, // inputs
        long oN[], long oB[], flt oA[], flt ob[], flt oc[], flt &ov, flt oc2[], flt &ov2                        // outputs
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
    ob[i] = max(b[i] - A[index(i, e, n)] * ob[l], 0);
    for (long j = 0; j < n; j++) {
      if (j == e) continue;
      oA[index(i, j, n)] = A[index(i, j, n)] - A[index(i, e, n)] * oA[index(l, j, n)];
    }
    oA[index(i, e, n)] = -A[index(i, e, n)] * oA[index(l, e, n)];
  }

  // Compute the objective function
  ov = v + c[e] * ob[l];
  if (c2 != nullptr) ov2 = v2 + c2[e] * ob[l];
  for (long j = 0; j < n; j++) {
    if (j == e) continue;
    oc[j] = c[j] - c[e] * oA[index(l, j, n)];
    if (c2 != nullptr) oc2[j] = c2[j] - c2[e] * oA[index(l, j, n)];
  }
  oc[e] = -c[e] * oA[index(l, e, n)];
  if (c2 != nullptr) oc2[e] = -c2[e] * oA[index(l,e,n)];

  // Compute new sets of basic and nonbasic variables
  if (oB != B || N != oN) {
    for (long i = 0; i < n; i++) {
      oB[i] = B[i];
    }
    for (long j = 0; j < n; j++) {
      oN[j] = N[j];
    }
  }
  swap(oB[l], oN[e]);
}

std::string simplex_slack(
        long N[], long B[], flt A[], flt b[], flt c[], flt v, long n, long m, flt maxZ, flt c2[], flt &v2, // inputs
        flt &z, long &iters                                                             // outputs
) {
  if (PRINT)
    print_slack(N, B, A, b, c, v, n, m);
  long e = 0;
  for (long j = 1; j < n; j++)
    e = c[j]>c[e] ? j : e;
  while (c[e] > EPS) { // While there exists j with c_j > 0 && v <= maxZ - EPS

    if (DEBUG && iters % 100 == 0) {
      printf("k=%li, c[e]=%.17g, v=%.17g, maxZ-EPS=%.17g, v2=%17.g\n", iters, c[e], v, maxZ - EPS, v2);
    }

    // Find i in B minimizing b_i / a_ie
    long l = -1;
    for (long i = 0; i < m; i++) {
      if (A[index(i, e, n)] > EPS &&
          (l == -1 || b[l] / A[index(l, e, n)] > b[i] / A[index(i, e, n)]))
        l = i;
    }
    if (l == -1) {
      char buffer [100];
      snprintf(buffer, 100, "unbounded in x%li (col %li)", N[e], e);
      z = INF;
      return std::string(buffer);
    }

    if (DEBUG) printf("x%li (col %li) enters; x%li (row %li) leaves. ce=%.17g\n", N[e], e, B[l], l, c[e]);

    pivot(N, B, A, b, c, v, c2, v2, l, e, n, m,        // inputs
          N, B, A, b, c, v, c2, v2);                   // outputs
    iters++;

    if (PRINT) print_slack(N, B, A, b, c, v, n, m);

    e = 0;
    for (long j = 1; j < n; j++)
      e = c[j]>c[e] ? j : e;
  }
  z = v;

  if (DEBUG) printf("k=%li, c[e]=%.17g, v=%.17g, maxZ-EPS=%.17g, v2=%17.g\n", iters, c[e], v, maxZ - EPS, v2);
  return "success";
}

std::string simplex(
        flt A[], flt b[], flt c[], long n, long m,
        flt &z, long* B, long &iters
) {
  long *N = new long[n];
  iters = 0;

  // PHASE 1: Find initial solution.

  // Find k minimizing b[i]
  long k = 0;
  for (long i = 1; i < m; i++)
    k = b[i] < b[k] ? i : k;
  if (b[k] >= 0) { // We can use the canonical slack form initially.
    for (long j = 0; j < n; j++)
      N[j] = j;
    for (long i = 0; i < m; i++)
      B[i] = n + i;
    z = 0;
  } else { // We have to solve an auxiliary problem, to find an initial solution.
    flt *auxA = new flt[m * (n + 1)];
    flt *auxc = new flt[n + 1];
    flt *auxb = new flt[m];
    long *auxN = new long[n + 1];
    long *auxB = new long[m];

    for (long j = 0; j < n; j++) {
      auxc[j] = 0.;
      auxN[j] = j;
    }
    auxc[n] = -1.;
    auxN[n] = n+m;

    c[n] = 0;

    for (long i = 0; i < m; i++) {
      auxB[i] = n + i;
      auxb[i] = b[i];
    }
    for (long i = 0; i < m; i++) {
      for (long j = 0; j < n; j++) {
        auxA[index(i, j, n + 1)] = A[index(i, j, n)];
      }
      auxA[index(i,n,n+1)] = -1.;
    }

    if (DEBUG) printf("Before auxilary pivot of e=n and l=%li\n", k);

    flt auxv = 0;
    if (PRINT) { printf("Auxilary problem:\n"); print_slack(auxN, auxB, auxA, auxb, auxc, auxv, n+1, m); }

    pivot(auxN, auxB, auxA, auxb, auxc, auxv, c, z, k, n, n + 1, m,
          auxN, auxB, auxA, auxb, auxc, auxv, c, z);
    iters++;


    if (DEBUG) {
      printf("After auxilary pivot of e=n and l=%li\n", k);
      for (long i = 0; i < m; i++) {
        if (auxb[i] < 0.) {
          printf("b%li = %lf\n", i, auxb[i]);
        }
      }
    }

    flt auxz;
    std::string result = simplex_slack(auxN, auxB, auxA, auxb, auxc, auxv, n + 1, m, 0., c, z,
                                       auxv, iters);
    auxz = auxv;
    if (result != "success") {
      printf("Could not solve the auxiliary problem for finding an initial value.\n"
             "This should never happen. Instead: %s\n", result.c_str());
      return "error";
    }
    if (auxz < -EPS) {
      if (DEBUG) printf("The problem was found to be infeasible.");
      return "infeasible";
    }


    if (DEBUG) {
      printf("Before degenerate pivot\n");
      for (long i = 0; i < m; i++) {
        if (auxb[i] < 0) {
          printf("b%li = %lf\n", i, auxb[i]);
        }
      }
    }

    // Check whether n+m is a basis variable. If so, do a non-degenerate pivot.
    for (long l = 0; l < m; l++) {
      if (auxB[l] == n + m){
        long e = 0;
        for (long j = 1; j < n; j++) {
          if (std::abs(auxA[index(l, j, n + 1)]) > std::abs(auxA[index(l, e, n + 1)]))
            e = j;
        }
        if (DEBUG) printf("Erasing the pivot to keep solution feasible: b[l]=%.17g\n", auxb[l]);
        auxb[l] = 0;
        pivot(auxN, auxB, auxA, auxb, auxc, auxv, c, z, l, e, n + 1, m,
              auxN, auxB, auxA, auxb, auxc, auxv, c, z);
        iters++;
        break;
      }
    }

    if (DEBUG) {
      printf("After degenerate pivot\n");
      for (long i = 0; i < m; i++) {
        if (auxb[i] < 0) {
          printf("b%li = %lf\n", i, auxb[i]);
        }
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
      c[j > remCol ? j - 1 : j] = c[j];
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
  }

  // PHASE 2: Find optimal solution.
  return simplex_slack(N, B, A, b, c, z, n, m, INF, nullptr, z,
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

  std::string result = simplex_slack(N, B, A, b, c, v, n, m, INF, nullptr, z,
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

void testFromFile(int argc, char ** argv) {
  std::string inputStr;
  if (argc < 3) {
    char* inputFile = new char[100];
    printf("Please enter problem name:\n");
    if (scanf("%s", inputFile) != 1) printf("Entered problem name invalid!");
    std::string inputStr = std::string(inputFile);
  } else {
    inputStr = std::string(argv[2]);
  }

  long m, n;
  std::ifstream fileshape = std::ifstream(inputStr + "-shape.csv");
  if (!fileshape) std::cout << "Could not open file: " + inputStr + "-shape.csv";
  std::string cell;
  std::getline(fileshape, cell);
  m = std::stol(cell);
  std::getline(fileshape, cell);
  n = std::stol(cell);
  fileshape.close();

  flt * A = new flt[m*n];
  std::ifstream fileA = std::ifstream (inputStr + "-A.csv");
  if (!fileA) std::cout << "Could not open file: " + inputStr + "-A.csv";
  for (long i = 0; i < m; i++) {
    std::string cell;
    for (long j = 0; j < n; j++) {
      if(!std::getline(fileA, cell, ',')) printf("Couldn't read A");
      A[index(i,j,n)] = std::stod(cell);
    }
  }
  fileA.close();

  flt * b = new flt[m];
  std::ifstream fileb = std::ifstream (inputStr + "-b.csv");
  if (!fileb) std::cout << "Could not open file: " + inputStr + "-b.csv";
  for (long i = 0; i < m; i++) {
    if(!std::getline(fileb, cell)) printf("Couldn't read b");
    b[i] = std::stod(cell);
  }
  fileb.close();

  flt * c = new flt[n + 1];
  std::ifstream filec = std::ifstream (inputStr + "-c.csv");
  if (!filec) std::cout << "Could not open file: " + inputStr + "-c.csv";
  for (long j = 0; j < n; j++) {
    std::getline(filec, cell);
    c[j] = std::stod(cell);
  }
  filec.close();

  flt z;
  long * B = new long[m];

  if (argc < 3) printf("Starting Linear Optimization...\n");
  long iters = 0;
  auto start = std::chrono::high_resolution_clock::now();

  std::string result = simplex(A, b, c, n, m, z, B, iters);

  auto finish = std::chrono::high_resolution_clock::now();
  double ns = std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();
  double s = ((double) ns) * 1e-9;
  if (argc < 3) {
    printf("Finished Linear Optimization with result %s (and optimal value %lf)\n", result.c_str(), z);
    printf("This took %lf seconds\n", s);
  }
  else printf(R"({"m": %li, "n": %li, "name": "%s", "t": %lf, "k": %li, "z": %lf},)" "\n", m, n, inputStr.c_str(), s, iters, z);

  delete[] A;
  delete[] B;
  delete[] b;
  delete[] c;
}

void testFromRand (int argc, char ** argv) {
  long n;
  if (argc != 3) {
    printf("Choose matrix size nxn!\n");
    if (scanf("%li", &n) != 1) printf("Invalid input!");
  } else n = std::stol(argv[2]);
  long m = n;

  flt z;
  long * B = new long[m];

  std::uniform_real_distribution<flt> unif(-1., 1.);
  std::default_random_engine re(12345);

  flt * A = new flt[m*n];
  flt * b = new flt[m];
  flt * c = new flt[n + 1];
  for (long i = 0; i < m; i++) {
    for (long j = 0; j < n; j++)
      A[index(i,j,n)] = unif(re);
    b[i] = (unif(re) + 1)/2;
  }
  for (long j = 0; j < n; j++)
    c[j] = unif(re);


  if (argc != 3) printf("Starting Linear Optimization...\n");
  long iters = 0;
  auto start = std::chrono::high_resolution_clock::now();

  std::string result = simplex(A, b, c, n, m, z, B, iters);

  auto finish = std::chrono::high_resolution_clock::now();
  if (argc != 3) printf("Finished Linear Optimization with result %s and optimal value %lf\n", result.c_str(), z);
  double ns = std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count();
  double s = ((double) ns) * 1e-9;
  if (argc != 3) printf("This took %lf seconds and %li iterations\n", s, iters);
  else printf(R"({"n": %li, "t": %lf, "k": %li, "z": %lf},)" "\n", n, s, iters, z);

  delete[] A;
  delete[] B;
  delete[] b;
  delete[] c;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    printf("Add file|rand as command line argument!\n");
    exit(EXIT_FAILURE);
  } else if (strcmp("file", argv[1]) == 0) {
    testFromFile(argc, argv);
  } else if (strcmp("rand", argv[1]) == 0) {
    testFromRand(argc, argv);
  }
  return 0;
}
