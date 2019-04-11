#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>
#include <assert.h>
#include "athread.h"
#include "swBLAS_def.h"

//#define DBL_EPSILON 1e-10
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

#define lda_ 64
#define ldc_ 64
#define ndrow1_ 16
#define ndrow3_ (lda-ndrow1)
#define wdt_ 56

extern SLAVE_FUN(fusion_dgemm_dsyrk_single_slave)();



void init_fusion(double* a, double* c, double* c2, int lda, int ldc, int wdt, int ndrow1){
    int i;
    srand((unsigned int) 1);
    for (i = 0; i < lda*wdt; i++){
        a[i] = i%100 * 0.01;
    }
    for (i = 0; i < ldc*ndrow1; i++){
        c[i] = 0.;
        c2[i] = 0.;
    }
    return;
}
void check_fusion(double* c, double* c2, int ldc, int m, int n){
    int i, j;
    int offset;
    int err_cnt = 0;
    int total_cnt = 0;
    for (j = 0; j < n; j++){
        for (i = 0; i < m; i++){
            offset = j*ldc + i;
            if( (c[offset] - c2[offset]) > DBL_EPSILON || 
                (c[offset] - c2[offset]) < -DBL_EPSILON ){
                err_cnt++;
                printf("(%d,%d) c=%lf c2=%lf\n", j,i,c[offset],c2[offset]);
            }
            total_cnt++;
        }
    }
    printf("%d error out of shape %d %d total %d\n", err_cnt, m, n, total_cnt);
    return;
}

int main(int argc, char** argv){
    athread_init();
    struct timeval t1, t2;
    double blasTime = 0., myTime = 0.;

    int lda, ldc, ndrow1, ndrow3, wdt;

    lda = lda_;
    ldc = ldc_;
    ndrow1 = ndrow1_;
    ndrow3 = ndrow3_;
    wdt = wdt_;

    assert(lda >= ndrow1+ndrow3);
    assert(ldc >= ndrow1+ndrow3);

    double* input_a = (double*)malloc(sizeof(double)*lda*wdt);
    double* input_c = (double*)malloc(sizeof(double)*ldc*ndrow1);
    double* input_c2 = (double*)malloc(sizeof(double)*ldc*ndrow1);

    char TRANSA, TRANSB;
    char UPLO, Tran;
    double one = 1., zero = 0.;
    init_fusion(input_a, input_c, input_c2, lda, ldc, wdt, ndrow1);
int ii;
for (ii = 0; ii < 10; ii++){
    UPLO = 'L';
    Tran = 'N';
    TRANSA = 'N';
    TRANSB = 'C';
    gettimeofday(&t1, NULL);
    dsyrk_(&UPLO,&Tran,&ndrow1,&wdt,&one,input_a,&lda,&zero,input_c,&ldc);
    dgemm_(&TRANSA, &TRANSB, &ndrow3, &ndrow1, &wdt, &one, \
            input_a+ndrow1, &lda, input_a, &lda, &zero, input_c+ndrow1, &ldc);
    gettimeofday(&t2, NULL);
    blasTime += TIME(t1, t2);

    fusionDgemmDsyrkParam param;
    param.lda = lda;
    param.ldc = ldc;
    param.ndrow1 = ndrow1;
    param.ndrow3 = ndrow3;
    param.wdt = wdt;
    param.a = input_a;
    param.c = input_c2;

    printf("BEGIN single_slave\n");

    gettimeofday(&t1, NULL);
    athread_spawn(fusion_dgemm_dsyrk_single_slave, &param);
    athread_join();
    gettimeofday(&t2, NULL);
    myTime += TIME(t1, t2);
}
    printf("blasTime is %lf \n", blasTime);
    printf("myTime is %lf \n", myTime);

    check_fusion(input_c, input_c2, ldc, (ndrow1+ndrow3), ndrow1);

    free(input_a);
    free(input_c);
    free(input_c2);

    return 0;
}


