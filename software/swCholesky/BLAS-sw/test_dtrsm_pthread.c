#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include "athread.h"
#include "swBLAS_def.h"

#include <pthread.h>
#include <unistd.h>

#define EPSILON 1e-6
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

extern SLAVE_FUN(dsyrk_single_slave)(void*);

/*
#define lda_ 1463
#define ldb_ 1463
#define m_ 1120
#define n_ 343
*/

/*
#define lda_ 48
#define ldb_ 48
#define m_ 35
#define n_ 13
*/


#define lda_ 6400
#define ldb_ 6400
#define m_ 6400
#define n_ 6400


/*
#define lda_ 64
#define ldb_ 64
#define m_ 64
#define n_ 64
*/

#define NUM_CG 4

void init_dtrsm(double* a, double* b, double* b2, int m, int n, int lda, int ldb){
    int i;
    srand((unsigned int) 1);
    for (i = 0; i < lda*n; i++){
        a[i] = i%100 * 0.01;
    }
    for (i = 0; i < n; i++){
        a[i*lda+i] += 2.;
        //printf("a[%d %d] = %lf\n", i,i,a[i*lda+i]);
    }
    for (i = 0; i < ldb*n; i++){
        b[i] = 1. + i%100*0.01;
        b2[i] = b[i];
    }
    return;
}
void check_dtrsm(double* c, double* c2, int ldc, int m, int n){
    int i, j;
    int offset;
    int err_cnt = 0;
    int total_cnt = 0;
    for (i = 0; i < n; i++){
        for (j = 0; j < m; j++){
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



int pthread_main(){
    struct timeval t1, t2;
    double blasTime, myTime;

    int lda, ldb, m, n;

    lda = lda_;
    ldb = ldb_;
    m = m_;
    n = n_;

    printf("lda=%d ldb=%d m=%d n=%d\n", lda,ldb,m,n);

    double* input_a = (double*)malloc(sizeof(double)*lda*n);
    double* input_b = (double*)malloc(sizeof(double)*ldb*n);
    double* input_b2 = (double*)malloc(sizeof(double)*ldb*n);

    double one = 1., zero = 0.;
    init_dtrsm(input_a, input_b, input_b2, m, n, lda, ldb);

    printf("@FINISH init\n");

    char SIDE = 'R';
    char UPLO = 'L';
    char Trana = 'C';
    char DIAG = 'N';
    gettimeofday(&t1, NULL);
    //dtrsm_(&SIDE, &UPLO, &Trana, &DIAG, &m, &n,&one, \
                        input_a,&lda,input_b,&ldb);
    gettimeofday(&t2, NULL);
    blasTime = TIME(t1, t2);

    printf("@FINISH xmath\n");

    dtrsmParam param;
    param.m = m;
    param.n = n;
    param.lda = lda;
    param.ldb = ldb;
    param.a = input_a;
    param.b = input_b2;

    printf("BEGIN single_slave\n");

    gettimeofday(&t1, NULL);
    athread_init();
    athread_init();
    athread_init();
    athread_init();
    athread_init();
    //athread_spawn(dsyrk_single_slave, &param);
    //athread_join();
    //athread_halt();
    //dtrsm_single(m, n, input_a, lda, input_b2, ldb);
    dtrsm_cluster_master(m, n, one, input_a, lda, input_b2, ldb);
    gettimeofday(&t2, NULL);
    myTime = TIME(t1, t2);

    printf("blasTime is %lf \n", blasTime);
    printf("myTime is %lf \n", myTime);
    printf("SPEEDUP is %lf \n", blasTime/myTime);

    //check_dtrsm(input_b, input_b2, ldb, m, n);

    printf("Begin FREE\n");
    free(input_a);
    printf("FREE a\n");
    free(input_b);
    printf("FREE b\n");
    free(input_b2);
    printf("FREE b2\n");

    return 0;
}

int main(){
    pthread_t pthread_handler[NUM_CG];

    int iii;
    for(iii = 0; iii < NUM_CG; iii++){

        int rc = pthread_create(&pthread_handler[iii], NULL, pthread_main, NULL);
        if (rc != 0)
            printf("ERROR ##FINISH pthread create%d, rc=%d\n", iii, rc);
        else
            printf("##FINISH pthread create%d, rc=%d\n", iii, rc);
    }
    for(iii = 0; iii < NUM_CG; iii++){
        pthread_join(pthread_handler[iii], NULL);
        printf("##FINISH pthread join%d\n", iii);
    }


    printf("##FINISH pthread\n");
    return 0;
}
