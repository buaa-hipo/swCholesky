#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <athread.h>
#include "cblas.h"
#include "swBLAS_def.h"
#include "swBLAS_fun.h"

//ERROR
void dpotrf2_cluster_master(const int n, double* a, const int lda, int* info){

    const double one = 1.0;
    const double mone = -1.0;
    const double zero = 0.;

    *info = 0;

    int iinfo;
    int n1, n2;

    if(n == 0) return;

    if (n == 1){
        if(a[0] <= 0. || isnan(a[0]) ){
            printf("n==1 %lf\n", a[0]);
            *info = 1;
            return;
        }

        a[0] = sqrt(a[0]);
    }
    else{
        n1 = n/2;
        n2 = n - n1;

        printf("A11 n = %d n1 = %d, n2 = %d\n", n, n1, n2);
        dpotrf2_cluster_master(n1, a, lda, &iinfo);
        if (iinfo != 0){
            *info = iinfo;
            return;
        }

        dtrsm_cluster_master(n2, n1, one, a, lda, a+(n1), lda);
        //cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, n2, n1, one, a, lda, a+(n1), lda);
        printf("A22 n = %d n1 = %d, n2 = %d\n", n, n1, n2);
        dsyrk_cluster_master(n2, n1, mone, a+(n1)*lda, lda, one, a+(n1)*lda+(n1), lda);
        dpotrf2_cluster_master(n2, a+(n1)*lda+(n1), lda, &iinfo);
        if(iinfo != 0){
            *info = iinfo + n1;
            return;
        }

    }

}
