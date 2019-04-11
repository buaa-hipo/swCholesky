#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <athread.h>
#include "swBLAS_def.h"
#include "swBLAS_fun.h"
#include "cblas.h"

#define DEBUG_VERBOSE
#undef DEBUG_VERBOSE

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


extern SLAVE_FUN(dtrsm_cluster_slave)();
extern SLAVE_FUN(dtrsm_pad_b)();
extern SLAVE_FUN(dtrsm_depad_b)();

void dtrsm_cluster_master(const int m, const int n, const double alpha, const double* a, const int lda, double* b, const int ldb){

    if(m == 0 || n == 0) return;

    assert(alpha == 1.0);

#ifndef PURE
    if (n <= 32 && m <= 64){
        dtrsm_single(m,n,alpha,a,lda,b,ldb); 
        return;
    }
#endif

#ifndef CHOLESKY_BLAS
    char SIDE = 'R';
    char UPLO = 'L';
    char Trana = 'C';
    char DIAG = 'N';
#endif

    // determine the block size
    int blkN = MIN(256, n);
    int numN = (n-1)/blkN+1;
    int blkM_perCPE = MIN( (56*1024)/(blkN*8)-2, (m-1)/64+1 );
    blkM_perCPE = ((blkM_perCPE-1)/4+1)*4;

    while (blkM_perCPE > 8){
        int ldm_use = sizeof(double)*( blkN*(2 + blkM_perCPE) )/1024;
        #ifdef DEBUG_VERBOSE
        if (ldm_use < 60) {
            printf("ldm_use=%d\n", ldm_use);
            break;
        }
        #endif
        blkM_perCPE -= 4;
    }
    int blkM = blkM_perCPE*64;
    int numM = (m-1)/blkM+1; //ceiling
    int Ms = m/blkM*blkM;
    int Me = (m+blkM-1)/blkM*blkM;

    double *bp_ = (double*)malloc(sizeof(double)*n*blkM + 32);
    double *bp = ALIGNED(bp_);
#ifdef DEBUG_VERBOSE
    printf("m=%d, Ms=%d, Me=%d, blkM=%d, numM=%d\n", m,Ms,Me,blkM,numM);
#endif
    dtrsmPadParam paramPad;
    paramPad.m = m;
    paramPad.n = n;
    paramPad.ldb = ldb;
    paramPad.blkM = blkM;
    paramPad.b = b;
    paramPad.bp = bp;

    if(Me > m){
#ifdef DEBUG_VERBOSE
        printf("padding\n");
#endif
        athread_spawn(dtrsm_pad_b, &paramPad);
        athread_join();
#ifdef DEBUG_VERBOSE
        printf("##END padding\n");
#endif
    }

    int i, j, cnt=0;

#ifdef DEBUG_VERBOSE
    if (Me > m){
    for (j = 0; j < n; j++){
        for (i = Ms; i < Me; i++){
            int off = j*blkM + i-Ms;
            if (i >= m){
                if (bp[off] != 0){
                    printf("bp[%d %d]=%lf, should=%lf\n", i,j,bp[off], 0);
                    cnt++;
                }
            }
            else {
                double val = 1.+ (j*ldb+i)%100*0.01 ;
                if (bp[off] != val || bp[off] != b[j*ldb+i]){
                    printf("bp[%d %d]=%lf, b=%lf, should=%lf\n", i,j,bp[off], b[j*ldb+i], val);
                    printf("%d*%d+%d = %d\n", j,ldb,i, j*ldb+i);
                    cnt++;
                }
            }
        }
    }
    printf("ERROR CNT = %d\n", cnt);
    }
    //exit(0);
#endif
    

    int cN, cM;

    const double one = 1.0;
    const double mone = -1.0;

    dtrsmParam param;
    for (cN = 0; cN < numN; cN++){
        int curblkN = MIN(blkN, n - cN*blkN);
        for (cM = 0; cM < numM; cM++){
            if((cM+1)*blkM <= m){
                param.m = blkM;
                param.n = curblkN;
                param.lda = lda;
                param.ldb = ldb;
                param.a = (void*)(a + cN*blkN*lda + cN*blkN);
                param.b = (void*)(b + cN*blkN*ldb + cM*blkM);
            }
            else{
                param.m = blkM;
                param.n = curblkN;
                param.lda = lda;
                param.ldb = blkM;
                param.a = (void*)(a + cN*blkN*lda + cN*blkN);
                param.b = (void*)(bp + cN*blkN*blkM);
            }

#ifdef DEBUG_VERBOSE
            printf("%d %d | %d %d | %d %d |\n", param.n, param.m, numN, numM, cN, cM);
            printf("a=%p, b=%p, param.a=%p, param.b=%p\n", a, b, param.a, param.b);
#endif

#ifdef CHOLESKY_BLAS
            athread_spawn(dtrsm_cluster_slave, &param);
            athread_join();
#else
            dtrsm_(&SIDE, &UPLO, &Trana, &DIAG, &param.m, &param.n,&one,
                        param.a,&param.lda,param.b,&param.ldb);
#endif

#ifdef DEBUG_VERBOSE
            printf("END\n");
#endif
        }

        if (n - cN*blkN - curblkN > 0){
#ifdef DEBUG_VERBOSE
            printf("NEED GEMM\n");
#endif
/*          
            if (Ms > 0){
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, \
                        Ms, (n-cN*blkN-curblkN), curblkN, mone, \
                        b + cN*blkN*ldb, ldb, \
                        a + cN*blkN*lda + cN*blkN + curblkN, lda, \
                        one, b + (cN*blkN+curblkN)*ldb, ldb);
            }
            if (Me > m){
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, \
                        blkM, (n-cN*blkN-curblkN), curblkN, mone, \
                        bp + cN*blkN*blkM, blkM, \
                        a + cN*blkN*lda + cN*blkN + curblkN, lda, \
                        one, bp + (cN*blkN+curblkN)*blkM, blkM);
            }
*/                        
            if (Ms > 0){
                #ifdef DEBUG_VERBOSE
                printf("GEMM 1\n");
                #endif
                dgemm_cluster_master( \
                        Ms, (n-cN*blkN-curblkN), curblkN, mone, \
                        b + cN*blkN*ldb, ldb, \
                        a + cN*blkN*lda + cN*blkN + curblkN, lda, \
                        one, b + (cN*blkN+curblkN)*ldb, ldb);
            }
            if (Me > m){
                #ifdef DEBUG_VERBOSE
                printf("GEMM 2\n");
                #endif
                dgemm_cluster_master( \
                        blkM, (n-cN*blkN-curblkN), curblkN, mone, \
                        bp + cN*blkN*blkM, blkM, \
                        a + cN*blkN*lda + cN*blkN + curblkN, lda, \
                        one, bp + (cN*blkN+curblkN)*blkM, blkM);
            }
            //cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, \
                        (n-cN*blkN-curblkN), m, curblkN, mone, \
                        a + cN*blkN*lda + cN*blkN + curblkN, lda, \
                        b + cN*blkN*ldb, ldb, \
                        one, b + (cN*blkN+curblkN)*ldb, ldb);
           // dgemm_cluster_master( \
                        m, (n-cN*blkN-curblkN), curblkN, mone, \
                        b + cN*blkN*ldb, ldb, \
                        a + cN*blkN*lda + cN*blkN + curblkN, lda, \
                        one, b + (cN*blkN+curblkN)*ldb, ldb);
        }
    }

    if(Me > m){
#ifdef DEBUG_VERBOSE
        printf("depadding\n");
#endif
        athread_spawn(dtrsm_depad_b, &paramPad);
        athread_join();
#ifdef DEBUG_VERBOSE
        printf("##END depadding\n");
#endif
    }

    free(bp_);

    return;
}
