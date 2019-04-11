#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <athread.h>
#include "swBLAS_def.h"
#include "swBLAS_fun.h"
#include "cblas.h"


void dsyrk_cluster_master(const int n, const int k, const double alpha, const double* a, const int lda, const double beta, double* c, const int ldc){

    if ((n == 0) || (((alpha == 0.) || (k == 0)) && (beta == 1.0))) return;

#ifndef PURE
    if ( (n+k <= 64) || !(n >= 32 && k >= 32) ){
        dsyrk_single(n,k,alpha,a,lda,beta,c,ldc);
        return;
    }
#endif

    assert( (alpha == -1.0 && beta == 1.0) || 
            (alpha ==  1.0 && beta == 0.0) );
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);
    //cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);
    //cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, \
                n, k, alpha, a, lda, beta, c, ldc);
    //sw_cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);
    dgemm_cluster_master( \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);

    return;
}

void dsyrk_cluster_master_pure(const int n, const int k, const double alpha, const double* a, const int lda, const double beta, double* c, const int ldc){

    if ((n == 0) || (((alpha == 0.) || (k == 0)) && (beta == 1.0))) return;

    assert( (alpha == -1.0 && beta == 1.0) || 
            (alpha ==  1.0 && beta == 0.0) );
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);
    //cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);
    //cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, \
                n, k, alpha, a, lda, beta, c, ldc);
    //sw_cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);
    dgemm_cluster_master( \
                n, n, k, alpha, \
                a, lda, a, lda, \
                beta, c, ldc);

    return;
}
