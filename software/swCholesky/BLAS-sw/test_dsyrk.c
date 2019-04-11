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
#define k_ 30
#define n_ 30


extern SLAVE_FUN(dsyrk_single_slave)();



void init_dsyrk(double* a, double* c, double* c2, int lda, int ldc, int k, int n){
    int i;
    srand((unsigned int) 1);
    for (i = 0; i < lda*k; i++){
        a[i] = i%100 * 0.01;
    }
    for (i = 0; i < ldc*n; i++){
        c[i] = 0.;
        c2[i] = 0.;
    }
    return;
}
void check_dsyrk(double* c, double* c2, int ldc, int n){
    int i, j;
    int offset;
    int err_cnt = 0;
    int total_cnt = 0;
    for (j = 0; j < n; j++){
        for (i = j; i < n; i++){
            offset = j*ldc + i;
            if( (c[offset] - c2[offset]) > DBL_EPSILON || 
                (c[offset] - c2[offset]) < -DBL_EPSILON ){
                err_cnt++;
                printf("(%d,%d) c=%lf c2=%lf\n", j,i,c[offset],c2[offset]);
            }
            total_cnt++;
        }
    }
    printf("%d error out of rank %d total %d\n", err_cnt, n, total_cnt);
    return;
}

int main(int argc, char** argv){
    athread_init();
    struct timeval t1, t2;
    double blasTime, myTime;

    int lda, ldc, k, n;

    lda = lda_;
    ldc = ldc_;
    k = k_;
    n = n_;

    double* input_a = (double*)malloc(sizeof(double)*lda*k);
    double* input_c = (double*)malloc(sizeof(double)*ldc*n);
    double* input_c2 = (double*)malloc(sizeof(double)*ldc*n);

    char UPLO, TRAN;
    double one = 1., zero = 0.;
    init_dsyrk(input_a, input_c, input_c2, k, lda, ldc, n);

    UPLO = 'L';
    TRAN = 'N';
    gettimeofday(&t1, NULL);
    dsyrk_(&UPLO, &TRAN, &n, &k, &one, input_a, &lda, &zero, input_c, &ldc);
    gettimeofday(&t2, NULL);
    blasTime = TIME(t1, t2);

    dsyrkParam param;
    param.k = k;
    param.lda = lda;
    param.ldc = ldc;
    param.n = n;
    param.a = input_a;
    param.c = input_c2;

    printf("BEGIN single_slave\n");

    gettimeofday(&t1, NULL);
    athread_spawn(dsyrk_single_slave, &param);
    athread_join();
    //dsyrk_single(k, lda, ldc, n, input_a, input_c2);
    //dsyrk_cluster_master(n, k, one, input_a, lda, zero, input_c2, ldc);
    gettimeofday(&t2, NULL);
    myTime = TIME(t1, t2);

    printf("blasTime is %lf \n", blasTime);
    printf("myTime is %lf \n", myTime);

    check_dsyrk(input_c, input_c2, ldc, n);

    free(input_a);
    free(input_c);
    free(input_c2);

    return 0;
}


