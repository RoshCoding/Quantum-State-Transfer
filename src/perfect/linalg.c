#include "linalg.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void printdv(double *x, int n) {
  for (int i = 0; i < n; i++) {
    printf("%12.6lf \n", x[i]);
  }
}

void printdvf(complex double *x, int n, FILE *fileout) {
  //
  // Prints an n long double vector to a file already
  // open with FILE *fileout.
  //
  int i;

  for (i = 0; i < n; i = i + 1) {
    fprintf(fileout, "%12.6lf", creal(x[i]));
  }
  fprintf(fileout, "\n");
}

void fprintdvp(double complex *qt, FILE *fout, int n) {
  //
  // Prints an n long double vector times its complex conjugate to a file
  // already open with FILE *fileout.
  //
  int i;

  for (i = 0; i < n; i = i + 1) {
    fprintf(fout, "%12.6lf", creal(qt[i] * conj(qt[i])));
  }
  fprintf(fout, "\n");
}

void zero_vect(double complex *x, int n) {
  for (int i = 0; i < n; i++) {
    x[i] = 0.0;
  }
}

void fprintdv(double *x, FILE *fout, int n) {
  for (int i = 0; i < n; i++) {
    fprintf(fout, "%12.6lf", x[i]);
  }
}

void printzv(double complex *x, int n) {
  for (int i = 0; i < n; i++) {
    printf("%6.3lf + %6.3lfi \n", creal(x[i]), cimag(x[i]));
  }
}

void printdm(double *A, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%12.6lf", A[i + m * j]);
    }
    printf("\n");
  }
}

void printzm(double complex *A, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      printf("%6.3lf + %6.3lfi", creal(A[i + m * j]), cimag(A[i + m * j]));
    }
    printf("\n");
  }
}

void set_seed(int seed) { srand(seed); }

double rand_unif_betw(double a, double b) {
  double r = (double)rand() / RAND_MAX;
  return a + (b - a) * r;
}

void dfillrand(double *A, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A[i + m * j] = rand_unif_betw(-1, 1);
    }
  }
}

void zfillrand(double complex *A, int m, int n) {
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      A[i + m * j] = rand_unif_betw(-1, 1) + I * rand_unif_betw(-1, 1);
    }
  }
}

void fmatrix_like_richard(double complex *A, int n, double complex diag, double complex offdiag) {
  // Fill diagonal
  for (int i = 0; i < n; i++) {
    A[i + n * i] = diag;
  }
  // Fill off-diagonal
  for (int i = 0; i < n - 1; i++) {
    A[i + n * (i + 1)] = offdiag;
    A[i + 1 + n * i] = offdiag;
  }
}

void fmatrix_proper(double complex *A, int n) {
    // Christandl et al. perfect state transfer Hamiltonian
    // Diagonal elements are zero, off-diagonal couplings:
    // J_j = (1/2) * sqrt(j * (n - j))
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i + n * j] = 0.0;
        }
    }
    for (int j = 1; j < n; j++) {
        double J = 0.5 * sqrt(j * (n - j));
        A[(j - 1) + n * j] = J;
        A[j + n * (j - 1)] = J;
    }
}