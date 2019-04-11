#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include "athread.h"
#include "swBLAS_def.h"

#define EPSILON 1e-6
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

#define m_ 449
#define n_ 108
#define k_ 130

#define lda_ 842
#define ldc_ 557

// dsyrk_ + dgemm_  =  dsyrk_ + dgemm_cluster
// dsyrk_  =  dsyrk_cluster
// dgemm_  =  dgemm_cluster
// dsyrk_ + dgemm_  =  dsyrk_cluster + dgemm_ 

//ERROR ::
// dsyrk_ + dgemm_  !=  dsyrk_cluster + dgemm_cluster 
// dsyrk_cluster + dgemm_  !=  dsyrk_cluster + dgemm_cluster      //ERROR isnan(c2) 
// dsyrk_cluster + dgemm_  !=  dsyrk_cluster + dgemm_             //ERROR isnan(c2) 
// dsyrk_cluster           !=  dsyrk_cluster                      //ERROR isnan(c2) 
// dgemm_cluster           !=  dgemm_cluster                      //ERROR != 

void init_syrk_gemm(double* a, double* a2, double* c, double* c2, int lda, int ldc, int k, int m, int n){
    int i;
    srand((unsigned int) 1);
    for (i = 0; i < lda*k; i++){
        a[i] = 1. + i%100 * 0.01;
        a2[i] = a[i];
    }
    for (i = 0; i < ldc*n; i++){
        c[i] = 1. + i%100 * 0.01;
        c2[i] = c[i];
    }
    return;
}
void check_syrk_gemm(double* c, double* c2, int ldc, int m, int n){
    int i, j;
    int offset;
    int err_cnt = 0;
    int total_cnt = 0;
    for (j = 0; j < n; j++){
        for (i = j; i < m; i++){
            offset = j*ldc + i;
            assert( !isnan(c[offset]) );
            assert( !isnan(c2[offset]) );
            if( (c[offset] - c2[offset]) > EPSILON || 
                (c[offset] - c2[offset]) < -EPSILON ){
                err_cnt++;
                printf("(%d,%d) c=%lf c2=%lf\n", i,j,c[offset],c2[offset]);
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

    int lda, ldc, k, m, n;

    lda = lda_;
    ldc = ldc_;
    k = k_;
    m = m_;
    n = n_;

    double* input_a = (double*)malloc(sizeof(double)*lda*k);
    double* input_a2 = (double*)malloc(sizeof(double)*lda*k);
    double* input_c = (double*)malloc(sizeof(double)*ldc*n);
    double* input_c2 = (double*)malloc(sizeof(double)*ldc*n);

    char TRANSA, TRANSB;
    char UPLO, TRAN;
    double one = 1., zero = 0.;
    init_syrk_gemm(input_a, input_a2, input_c, input_c2, lda, ldc, k, m, n);
int ii;
for (ii = 0; ii < 1; ii++){
    UPLO = 'L';
    TRAN = 'N';
    TRANSA = 'N';
    TRANSB = 'C';
    gettimeofday(&t1, NULL);
    dsyrk_(&UPLO, &TRAN, &n, &k, &one, input_a, &lda, &zero, input_c, &ldc);
    //dsyrk_cluster_master(n, k, one, input_a, lda, zero, input_c, ldc);
    dgemm_(&TRANSA, &TRANSB, &m, &n, &k, &one, input_a+n, &lda, input_a, &lda, &zero, input_c+n, &ldc);
    //dgemm_cluster_master(m, n, k, one, input_a+n, lda, input_a, lda, zero, input_c+n, ldc);
    gettimeofday(&t2, NULL);
    blasTime += TIME(t1, t2);

    printf("BEGIN single_slave\n");

    gettimeofday(&t1, NULL);
    //dsyrk_(&UPLO, &TRAN, &n, &k, &one, input_a, &lda, &zero, input_c2, &ldc);
    dsyrk_cluster_master(n, k, one, input_a, lda, zero, input_c2, ldc);
    //dgemm_(&TRANSA, &TRANSB, &m, &n, &k, &one, input_a+n, &lda, input_a, &lda, &zero, input_c2+n, &ldc);
    dgemm_cluster_master(m, n, k, one, input_a+n, lda, input_a, lda, zero, input_c2+n, ldc);
    gettimeofday(&t2, NULL);
    myTime += TIME(t1, t2);
}
    printf("blasTime is %lf \n", blasTime);
    printf("myTime is %lf \n", myTime);

    check_syrk_gemm(input_c, input_c2, ldc, m+n, n);

    free(input_a);
    free(input_a2);
    free(input_c);
    free(input_c2);

    return 0;
}


