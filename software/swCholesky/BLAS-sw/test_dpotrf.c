#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include "athread.h"
#include "swBLAS_def.h"

#define EPSILON 1e-6
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

#define lda_ 1463
#define n_ 343





void init_dpotrf(double* a, double* a2, int n, int lda){
    int i;
    srand((unsigned int) 1);
    for (i = 0; i < lda*n; i++){
        a[i]  = 0.2 + i%100*0.001;
        a2[i] = a[i];
    }
    for (i = 0; i < n; i++){
        a[i*lda+i] += 20.;
        a2[i*lda+i] = a[i*lda+i];
        //printf("a[%d %d] = %lf\n", i,i,a[i*lda+i]);
    }
    return;
}
void check_dpotrf(double* c, double* c2, int ldc, int m, int n){
    int i, j;
    int offset;
    int err_cnt = 0;
    int total_cnt = 0;
    for (i = 0; i < n; i++){
        for (j = i; j < m; j++){
            offset = i*ldc + j;
            if( fabs(c[offset] - c2[offset]) > EPSILON ){ 
                err_cnt++;
                if(err_cnt < 100)printf("(%d,%d) c=%lf c2=%lf\n", j,i,c[offset],c2[offset]);
            }
            total_cnt++;
        }
    }
    printf("%d error out of rank %d %d total %d\n", err_cnt, m, n, total_cnt);
    return;
}

int main(int argc, char** argv){
    athread_init();
    struct timeval t1, t2;
    double blasTime, myTime;

    int lda, ldb, m, n;

    lda = lda_;
    n = n_;

    printf("lda=%d n=%d\n", lda,n);

    double* input_a = (double*)malloc(sizeof(double)*lda*n);
    double* input_a2 = (double*)malloc(sizeof(double)*lda*n);

    double one = 1., zero = 0.;
    init_dpotrf(input_a, input_a2, n, lda);

    printf("@FINISH init\n");

    char UPLO = 'L';
    int info1 = 0;
    gettimeofday(&t1, NULL);
    dpotrf_(&UPLO, &n, input_a, &lda, &info1);
    gettimeofday(&t2, NULL);
    blasTime = TIME(t1, t2);

    printf("@FINISH xmath\n");

    printf("BEGIN single_slave\n");

    int info2 = 0;
    gettimeofday(&t1, NULL);
    dpotrf_cluster_master(n, input_a2, lda, &info2);
    gettimeofday(&t2, NULL);
    myTime = TIME(t1, t2);

    printf("blasTime is %lf , info=%d\n", blasTime, info1);
    printf("myTime is %lf , info=%d\n", myTime, info2);
    printf("SPEEDUP is %lf \n", blasTime/myTime);

    check_dpotrf(input_a, input_a2, lda, n, n);

    printf("Begin FREE\n");
    free(input_a);
    printf("FREE a\n");
    free(input_a2);
    printf("FREE a2\n");

    return 0;
}


