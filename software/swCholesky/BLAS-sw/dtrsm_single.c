#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include "swBLAS_def.h"

/*
      SUBROUTINE DTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),B(LDB,*)
*/

//      SIDE = 'R'; UPLO = 'L'; Trana = 'C'; DIAG = 'N'
//      ALPHA = 1

//      a(n*n) lda >= n     b(m*n) ldb >= m


void dtrsm_single(const int m, const int n, const double alpha, const double* a, const int lda, double* b, const int ldb){
    if (m == 0 || n == 0) return;
    assert(alpha == 1.0);
    assert(lda >= n && ldb >= m);

    int k, j, i;
    double temp;


    for (k = 0; k < n; k++){
        //printf("[k = %d]\n", k);
        //assert(a[BLAS_OFFSET(k,k,lda)] != 0.);
        temp = 1./a[BLAS_OFFSET(k,k,lda)];
        for (i = 0; i < m; i++){
            b[BLAS_OFFSET(i,k,ldb)] *= temp;
        }
        //printf("\t 1\n");

        for (j = k+1; j < n; j++){
            //printf("[k = %d] [j = %d]\n", k, j);
            if(a[BLAS_OFFSET(j,k,lda)] != 0.){
                temp = a[BLAS_OFFSET(j,k,lda)];
                for (i = 0; i < m; i++){
                    b[BLAS_OFFSET(i,j,ldb)] -= temp*b[BLAS_OFFSET(i,k,ldb)];
                }
            }
        }
    }
    return;
}
