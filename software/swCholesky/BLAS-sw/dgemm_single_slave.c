#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <slave.h>
#include <dma.h>
#include "swBLAS_def.h"

void dgemm_single_slave_old(void* ptr){
    const int my_id = athread_get_id(-1);
    if (my_id != 0) return;
    volatile unsigned int get_reply = 0, put_reply = 0;

    dgemmParam* param = (dgemmParam*)ptr;

    int lda = param->lda;
    int ldb = param->ldb;
    int ldc = param->ldc;
    int k = param->k;
    int m = param->m;
    int n = param->n;
    double* a = param->a;
    double* b = param->b;
    double* c = param->c;

    //printf("BEGIN\n");
    //printf("%d %d %d %d %d %d\n", lda, ldb, ldc, k, m, n);

    int j, l, i;
    double temp;

    double buffer[MAX_SIZE_DOUBLE];
    double* c_local = buffer;
    double* a_local = c_local + m;
    double* b_local = a_local + k*m;

    assert(m+k*m+k*n <= MAX_SIZE_DOUBLE && "cache size");

    for (l = 0; l < k; l++){
        get_reply = 0;
        athread_get(PE_MODE, a+l*lda, a_local + l*m, 
                    sizeof(double)*m, &get_reply, 0,0,0);
        athread_get(PE_MODE, b+l*ldb, b_local + l*n, 
                    sizeof(double)*n, &get_reply, 0,0,0);
        while(get_reply != 2);
        //dma_wait(&get_reply, 2);
    }

    for (j = 0; j < n; j++){
        for (i = 0; i < m; i++){
            temp = 0.;
            for (l = 0; l < k; l++){
                temp += a_local[l*m + i] * b_local[l*n + j];
            }
            c_local[i] = temp;
        }
        put_reply = 0;
        athread_put(PE_MODE, c_local, c+j*ldc,
                    sizeof(double)*m, &put_reply, 0,0);
        while(put_reply != 1);
        //dma_wait(&put_reply, 1);
    }      

    
    return;
}


void dgemm_single_slave(void* ptr){
    const int my_id = athread_get_id(-1);
    if (my_id != 0) return;
    volatile unsigned int get_reply = 0, put_reply = 0;

    dgemmParam* param = (dgemmParam*)ptr;

    int lda = param->lda;
    int ldb = param->ldb;
    int ldc = param->ldc;
    int k = param->k;
    int m = param->m;
    int n = param->n;
    double* a = param->a;
    double* b = param->b;
    double* c = param->c;

    //printf("BEGIN\n");
    //printf("%d %d %d %d %d %d\n", lda, ldb, ldc, k, m, n);

    int j, l, i;
    double temp;

    double buffer[MAX_SIZE_DOUBLE];
    double* c_local = buffer;
    double* a_local = c_local + m;
    double* b_local = a_local + k*m;

    assert(m+k*m+k*n <= MAX_SIZE_DOUBLE && "cache size");

    volatile unsigned int replygetA=0, replygetB=0, replyputC=0;
    dma_desc dmagetA, dmagetB, dmaputC;
    dma_set_op(&dmagetA, DMA_GET);
    dma_set_mode(&dmagetA, PE_MODE);
    dma_set_reply(&dmagetA, &replygetA);
    dma_set_op(&dmagetB, DMA_GET);
    dma_set_mode(&dmagetB, PE_MODE);
    dma_set_reply(&dmagetB, &replygetB);
    dma_set_op(&dmaputC, DMA_PUT);
    dma_set_mode(&dmaputC, PE_MODE);
    dma_set_reply(&dmaputC, &replyputC);

    dma_set_size(&dmagetA, sizeof(double)*m*k);
    dma_set_bsize(&dmagetA, sizeof(double)*m);
    dma_set_stepsize(&dmagetA, sizeof(double)*(lda-m));  

    dma_set_size(&dmagetB, sizeof(double)*n*k);
    dma_set_bsize(&dmagetB, sizeof(double)*n);
    dma_set_stepsize(&dmagetB, sizeof(double)*(ldb-n));  

    dma(dmagetA, (long)a, (long)a_local);
    dma(dmagetB, (long)b, (long)b_local);
    dma_wait(&replygetA, 1);
    replygetA = 0;
    dma_wait(&replygetB, 1);
    replygetB = 0;

    dma_set_size(&dmaputC, sizeof(double)*m);
    dma_set_bsize(&dmaputC, 0);
    dma_set_stepsize(&dmaputC, 0);    


    for (j = 0; j < n; j++){
        for (i = 0; i < m; i++){
            temp = 0.;
            for (l = 0; l < k; l++){
                temp += a_local[l*m + i] * b_local[l*n + j];
            }
            c_local[i] = temp;
        }
        dma(dmaputC, (long)(c+j*ldc), (long)c_local);
        dma_wait(&replyputC, 1);
        replyputC = 0;
    }      

    
    return;
}