//
// zheev.c
// (Zomplex HErmitian EigenValues)
// Compile like:
// gcc zheev.c linalg.c -lm -llapack
//
#include "linalg.h"

void zheev_(char *jobz, char *uplo, int *n, double complex *a, int *lda,
            double *w, double complex *work, int *lwork, double *rwork,
            int *info);

void zgemm_(char *transa, char *transb, int *m, int *n, int *k,
            complex double *alpha, complex double *a, int *lda,
            complex double *b, int *ldb, complex double *beta,
            complex double *c, int *ldc);

void zgemv_(char *trans, int *m, int *n, complex double *alpha,
            complex double *a, int *lda, complex double *x, int *incx,
            complex double *beta, complex double *y, int *incy);

void timepass(int n, double complex *V, double *evals, double complex *qt,
              double complex *work1, double complex *work2, double t);

int main(void) {
  char jobz, uplo;
  char transa, transb, transc;
  int n, lda, lwork, info;
  double complex *A, *A_evec, *W, *work;
  double *w, *rwork;
  double alpha = 1;
  double beta = 1;
  double incx = 1;
  double incy = 1;

  // ------------ INITIALIZE (HERMITIAN) MATRIX A HERE ------------

  FILE *fileout;
  fileout = fopen("qst-out.txt", "w");

  n = 64; // dimension of A
  A = calloc(n * n, sizeof(double complex));

  fmatrix_like_richard(A, n, -13.6, -13.6, 10);

  // --------------------------------------------------------------

  jobz = 'V'; // compute both eigenvalues and eigenvectors

  uplo = 'L'; // lower triangular part of A is stored
              // 'U' for upper triangular part
              // it is okay to have both set, but you only
              // *need* just one

  lda = n; // matches dimension of A in 99% of cases
           // (rare when you need to set it to something else)

  // vector to hold the n (real) eigenvalues of A
  w = calloc(n, sizeof(double));

  lwork = (128 + 1) * n; // length of workspace array work
                         // (at least 2*n-1)

  // workspace array
  work = calloc(lwork, sizeof(double complex));

  // another workspace array (just size 3*n-2)
  rwork = calloc(3 * n - 2, sizeof(double));

  transa = 'A';
  transb = 'A';
  transc = 'C';

  printf("\nDiagonalizing the matrix A =\n");
  printzm(A, n, n);

  zheev_(&jobz, &uplo, &n, A, &lda, w, work, &lwork, rwork, &info);
  // Copy eigenvectors (columns of A) to A_evec
  A_evec = calloc(n * n, sizeof(double complex));
  for (int i = 0; i < n * n; i++) {
    A_evec[i] = A[i];
  }

  double complex *qt;
  qt = calloc(n, sizeof(complex double));
  qt[0] = 1;

  double complex *work1, *work2;
  work1 = calloc(2 * n, sizeof(complex double));
  work2 = calloc(2 * n, sizeof(complex double));
  double t = 0;
  printzv(qt, n);
  // timepass(n, A, W, qt, work1, work2, 0.1, p);
  // printdvf(qt, n, fileout);

  for (double t = 0; t < 50; t += 0.1) {
    timepass(n, A_evec, w, qt, work1, work2, t);
    fprintdvp(qt, fileout, n);
    zero_vect(qt, n);
    qt[0] = 1;
  }

  /*
  printf("\nAfter diagonalization:\n");
  printf("eigenvector[j] = column j of A\n");
  printf("eigenvalue [j] = w[j]\n\n");
  printf("A =\n");
  printzm(A, n, n);
  printf("\nqt =\n");
  printdv(qt, n);
  fprintdv(qt, fileout, n);
  printf("\n");
  */
  free(A);
  free(A_evec);
  free(w);
  free(qt);
  free(work);
  free(work1);
  free(work2);
  free(rwork);
  fclose(fileout);
  return 0;
}

void timepass(int n, double complex *V, double *evals, double complex *qt,
              double complex *work1, double complex *work2, double t) {
  // Implements qt = V * exp(-i E t) * V^H * qt
  char transn = 'N';
  char transc = 'C';
  double complex alpha = 1.0, beta = 0.0;
  int incx = 1, incy = 1;

  // work1 = V^H * qt
  zgemv_(&transc, &n, &n, &alpha, V, &n, qt, &incx, &beta, work1, &incy);

  // work2 = exp(-i E t) * work1
  for (int i = 0; i < n; i++) {
    work2[i] = cexp(-I * evals[i] * t) * work1[i];
  }

  // qt = V * work2
  zgemv_(&transn, &n, &n, &alpha, V, &n, work2, &incx, &beta, qt, &incy);
}
