#include <stdlib.h>
#include <stdio.h>
#include <float.h>

/*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
*     ..
*/

// TRANSA = 'N', TRANSB = 'C', ALPHA = 1, BETA = 0

//*           Form  C := alpha*A*B**T + beta*C
/*
void dgemm_single(  int k, int lda, int ldb, int ldc, int m, int n,
                    double* a, double* b, double* c){

    int j, i, l;
    double temp;
    for (j = 0; j < n; j++){
        for (i = 0; i < m; i++){
            temp = 0.;
            for (l = 0; l < k; l++){
                temp += a[l*lda + i] * b[l*ldb + j];
            }
            c[j*ldc + i] = temp;
        }
    }

    return;
}
*/

void dgemm_single(const int m, const int n, const int k, const double alpha, const double* a, const int lda, const double* b, const int ldb, const double beta, double* c, const int ldc){
    int j, i, l;
    double temp;

    if (m <= 0 || n <= 0 || k <= 0) return;

    for (j = 0; j < n; j++){
        if (beta == 0.){
            for (i = 0; i < m; i++){
                c[j*ldc + i] = 0.;
            }
        }
        else if (beta != 1.){
            for (i = 0; i < m; i++){
                c[j*ldc + i] *= beta;
            }
        }

        for (l = 0; l < k; l++){
            temp = alpha * b[l*ldb + j];
            for (i = 0; i < m; i++){
                c[j*ldc + i] += a[l*lda + i] * temp;
            }
        }
    }

    return;    
}