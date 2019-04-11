#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>

/*
 *     .. Scalar Arguments ..
       DOUBLE PRECISION alpha,beta
       INTEGER k,lda,ldc,n
       CHARACTER trans,uplo
 *     ..
 *     .. Array Arguments ..
       DOUBLE PRECISION a(lda,*),c(ldc,*)
*/
//  TRANS = 'N', UPLO = 'L', ALPHA = 1, BETA = 0
//void dsyrk_single(
//    int k, int lda, int ldc, int n,
//    double* a, double* c ){

void dsyrk_single(const int n, const int k, const double alpha, const double* a, const int lda, const double beta, double* c, const int ldc){
    
    if ((n == 0) || (((alpha == 0.) || (k == 0)) && (beta == 1.0))) return;
    
    assert(alpha == 1 && beta == 0);
    int j, l, i;
    double temp;

    for (j = 0; j < n; j++){
        //if (beta == 0.){
            for (i = j; i < n; i++){
                c[j*ldc + i] = 0.;
            }
        //}
        //else if (beta != 1.0){
        //    for (i = j; i < n; i++){
        //        c[j*ldc + i] = beta*c[j*ldc + i];
        //    }
        //}
        for (l = 0; l < k; l++){
            temp = a[l*lda + j];                
            if (temp != 0.){
                for (i = j; i < n; i++){
                    c[j*ldc + i] += temp*a[l*lda + i];
                }
            } 
        }
    }

    return;
}