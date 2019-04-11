//
// Created by kazem on 7/25/17.
//

#ifndef CHOLOPENMP_MYBLAS_H
#define CHOLOPENMP_MYBLAS_H

//#include <cmath>
#include <math.h>

void Cholesky_col(int n, int dim, double* a){
 double tmp = 0;
 int i,j,k;
 for (j = 0; j < dim; ++j) {
  for (k = 0; k < j; ++k) {
   tmp = a[k*n+j];
   for (i = j; i < dim; ++i) {
    a[j*n+i] = a[j*n+i] - a[k*n+i]*tmp;
   }
  }
  tmp = sqrt(a[j*n+j]);
  for (k = j+1; k < dim; ++k) {
   a[j*n+k] = a[j*n+k] / tmp;
  }
  a[j*n+j] = tmp;
 }
}

int lSolve_dense_col(int colSize,int col, double *M, double *rhs){
 int i,j;
 for (i = 0; i < col; ++i) {
  rhs[i*colSize]/=M[i*colSize+i];
  for (j = i+1; j < col; ++j) {
   rhs[j*colSize]-=M[i*colSize+j]*rhs[i*colSize];
  }
 }
 return 1;
}


#endif //CHOLOPENMP_MYBLAS_H
