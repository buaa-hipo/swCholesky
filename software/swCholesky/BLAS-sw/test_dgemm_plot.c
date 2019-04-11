#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>
#include <assert.h>
#include "athread.h"
#include "swBLAS_def.h"

//#define DBL_EPSILON 1e-10
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

#define lda_ 32
#define ldb_ 32
#define ldc_ 32
#define k_ 32
#define m_ 32
#define n_ 32


extern SLAVE_FUN(dgemm_single_slave)();
extern SLAVE_FUN(nofunc)();



void init_dgemm(double* a, double* b, double* c, double* c2, double* c3, int lda, int ldb, int ldc, int k, int m, int n){
    int i;
    srand((unsigned int) 1);
    for (i = 0; i < lda*k; i++){
        a[i] = i%100 * 0.01;
    }
    for (i = 0; i < lda*k; i++){
        b[i] = 1 + i%100 * 0.01;
    }
    for (i = 0; i < ldc*n; i++){
        c[i] = 0.;
        c2[i] = 0.;
        c3[i] = 0.;
    }
    return;
}
void check_dgemm(double* c, double* c2, int ldc, int m, int n){
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
    double blasTime = 0., mpeTime = 0., cpeTime = 0.;

    int lda, ldb, ldc, k, m, n;
/*
    lda = lda_;
    ldb = ldb_;
    ldc = ldc_;
    k = k_;
    m = m_;
    n = n_;
*/
int sz;
for(sz=2; sz <= 64; sz+=2){
    printf("#################size=%d\n", sz);
    lda=sz;ldb=sz;ldc=sz;
    k=sz;m=sz;n=sz;
    blasTime=0.;
    mpeTime=0.;
    cpeTime=0.;

    double* input_a = (double*)malloc(sizeof(double)*lda*k);
    double* input_b = (double*)malloc(sizeof(double)*ldb*k);
    double* input_c = (double*)malloc(sizeof(double)*ldc*n);
    double* input_c2 = (double*)malloc(sizeof(double)*ldc*n);
    double* input_c3 = (double*)malloc(sizeof(double)*ldc*n);

    char TRANSA, TRANSB;
    double one = 1., zero = 0.;
    init_dgemm(input_a, input_b, input_c, input_c2, input_c3, lda, ldb, ldc, k, m, n);
int ii;
int steps = 10;
for (ii = 0; ii < steps; ii++){
    TRANSA = 'N';
    TRANSB = 'C';
    gettimeofday(&t1, NULL);
    //dgemm_(&TRANSA, &TRANSB, &m, &n, &k, &one, input_a, &lda, input_b, &ldb, &zero, input_c, &ldc);
    dgemm_cluster_master_pure(m, n, k, one, input_a, lda, input_b, ldb, zero, input_c3, ldc);
    gettimeofday(&t2, NULL);
    blasTime += TIME(t1, t2);

    dgemmParam param;
    param.k = k;
    param.lda = lda;
    param.ldb = ldb;
    param.ldc = ldc;
    param.m = m;
    param.n = n;
    param.a = input_a;
    param.b = input_b;
    param.c = input_c2;

    //printf("BEGIN single_slave\n");
    gettimeofday(&t1, NULL);
    athread_spawn(dgemm_single_slave, &param);
    //athread_spawn(nofunc, 0);
    athread_join();
    gettimeofday(&t2, NULL);
    cpeTime += TIME(t1, t2);

    gettimeofday(&t1, NULL);
    dgemm_single(m, n, k, one, input_a, lda, input_b, ldb, zero, input_c3, ldc);
    gettimeofday(&t2, NULL);
    mpeTime += TIME(t1, t2);
}
    printf("blasTime is %lf \n", blasTime);
    printf("mpeTime is %lf \n", mpeTime);
    printf("cpeTime is %lf \n", cpeTime);

    //check_dgemm(input_c, input_c2, ldc, m, n);

    free(input_a);
    free(input_b);
    free(input_c);
    free(input_c2);
    free(input_c3);
}
    return 0;
}


