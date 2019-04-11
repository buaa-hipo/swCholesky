#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <athread.h>
#include "swBLAS_def.h"
#include "swBLAS_fun.h"
#include "cblas.h"


void dpotrf_cluster_master(const int n, double* a, const int lda, int* info){

    const int nb = 128;
    const double one = 1.0;
    const double mone = -1.0;
    const double zero = 0.;

    *info = 0;

    int j;
    char UPLO = 'L';

    if (n == 0) return;

    if (nb == 1 || nb >= n){
        //dpotrf_(&UPLO, &n, a, &lda, info);
        //dpotrf2_cluster_master(n, a, lda, info);
        dpotrf2_single(n, a, lda, info);
    }
    else{
        for (j = 0; j < n; j+=nb){
            //printf("j = %d\n", j);
            int jb = MIN(nb, n-j);
            //printf("\t dsyrk\n");
            dsyrk_cluster_master(jb, j, mone, a+j, lda, one, a+j*lda+j, lda);
            //cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, \
                        jb, j, mone, a+j, lda, one, a+j*lda+j, lda);
            //dpotrf_(&UPLO, &jb, a+j*lda+j, &lda, info);
            //dpotrf2_cluster_master(jb, a+j*lda+j, lda, info);
            //printf("\t dpotrf2\n");
            dpotrf2_single(jb, a+j*lda+j, lda, info);
            //dpotrf_(&UPLO, &jb, a+j*lda+j, &lda, info);
            if (*info != 0){
                *info = *info + j;
                printf("\t ERROR POTRF2, info=%d\n", *info);
                return;
            }

            if(j+jb < n){
                //printf("\t dgemm\n");
                dgemm_cluster_master(n-j-jb, jb, j, mone,
                                    a+(j+jb), lda, a+j, lda,
                                    one, a+j*lda+(j+jb), lda);
                //printf("\t dtrsm\n");
                dtrsm_cluster_master(n-j-jb, jb, one, a+j*lda+j, lda, a+j*lda+(j+jb), lda);
            }
            //printf("FINISH j = %d\n", j);
        }

    }
    return;

}
