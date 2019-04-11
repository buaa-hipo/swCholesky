#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <athread.h>
#include "swBLAS_def.h"
#include "swBLAS_fun.h"
#include "swblas.h"
#include "cblas.h"

#define DEBUG_VERBOSE
#undef DEBUG_VERBOSE

void dgemm_cluster_master(const int m, const int n, const int k, const double alpha, const double* a, const int lda, const double* b, const int ldb, const double beta, double* c, const int ldc){

    if (m <= 0 || n <= 0 || k <= 0) return;

#ifndef PURE
    //if ( (m+n+k < 64) || !(m >= 64 && n >= 32 && k >= 32) ){
    if ( (m < 64 && n < 32 && k < 32) ){
        dgemm_single(m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
        return;
    }
#endif

    assert(alpha == 1.0 || alpha == -1.0);
    assert(beta == 0. || beta == 1.0);
#ifndef CHOLESKY_BLAS
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, \
                m, n, k, alpha, \
                a, lda, b, ldb, \
                beta, c, ldc);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, m, k, alpha, \
                b, ldb, a, lda, \
                beta, c, ldc);
#else
    #ifdef DEBUG_VERBOSE
    printf("###swGEMM :: START\n");
    printf("m=%d, n=%d, k=%d\n", n,m,k);
    printf("lda=%d, ldb=%d, ldc=%d\n", ldb,lda,ldc);
    printf("alpha=%lf, beta=%lf\n", alpha,beta);
    #endif
    sw_cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, m, k, alpha, \
                b, ldb, a, lda, \
                beta, c, ldc);
    #ifdef DEBUG_VERBOSE
    printf("@@@swGEMM :: FINISH\n");
    #endif
#endif
    return;
}


void dgemm_cluster_master_pure(const int m, const int n, const int k, const double alpha, const double* a, const int lda, const double* b, const int ldb, const double beta, double* c, const int ldc){

    if (m <= 0 || n <= 0 || k <= 0) return;

    assert(alpha == 1.0 || alpha == -1.0);
    assert(beta == 0. || beta == 1.0);
#ifndef CHOLESKY_BLAS
    //cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, \
                m, n, k, alpha, \
                a, lda, b, ldb, \
                beta, c, ldc);
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, m, k, alpha, \
                b, ldb, a, lda, \
                beta, c, ldc);
#else
    #ifdef DEBUG_VERBOSE
    printf("###swGEMM :: START\n");
    printf("m=%d, n=%d, k=%d\n", n,m,k);
    printf("lda=%d, ldb=%d, ldc=%d\n", ldb,lda,ldc);
    printf("alpha=%lf, beta=%lf\n", alpha,beta);
    #endif
    sw_cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                n, m, k, alpha, \
                b, ldb, a, lda, \
                beta, c, ldc);
    #ifdef DEBUG_VERBOSE
    printf("@@@swGEMM :: FINISH\n");
    #endif
#endif
    return;
}
