#ifndef LINALG_H
#define LINALG_H

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void printdv(double *x, int n);
void printdvf(double complex *x, int n, FILE *fileout);
void fprintdvp(double complex *qt, FILE *fout, int n);
void zero_vect(double complex *x, int n);
void fprintdv(double *x, FILE *fout, int n);
void printzv(double complex *x, int n);
void printdm(double *A, int m, int n);
void printzm(double complex *A, int m, int n);
void set_seed(int seed);
double rand_unif_betw(double a, double b);
void dfillrand(double *A, int m, int n);
void zfillrand(double complex *A, int m, int n);
void fmatrix_like_richard(double complex *A, int n, double complex diag, double complex offdiag);
void fmatrix_proper(double complex *A, int n);

#endif