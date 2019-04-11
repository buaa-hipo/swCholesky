#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <slave.h>
#include <simd.h>
#include <dma.h>
#include <assert.h>
#include "swBLAS_def.h"

void dsyrk_single_slave_old(void* ptr){
    const int my_id = athread_get_id(-1);
    //if (my_id != 0) return;
    volatile unsigned int get_reply = 0, put_reply = 0;

    dsyrkParam* param = (dsyrkParam*)ptr;
    //get_reply = 0;
    //athread_get(PE_MODE, ptr, &param, 
    //            sizeof(dsyrkParam), &get_reply, 0,0,0);
    //while(get_reply != 1);

    // 在这里存在一些问题
    // athread_get param参数总是在while结束后 还没有获取完毕 
    // 当使用非-g模式时
    // 在这里借鉴了swdnn的方法，使用gload方式获取参数

    int j, l, i;
    double temp;

    //printf("%d %d %d %d\n\n", param.k, param.lda, param.ldc, param.n);
    int k = param->k;
    int lda = param->lda;
    int ldc = param->ldc;
    int n = param->n;
    double *a = param->a;
    double *c = param->c;

    //printf("%d %d %d %d\n", param.k, param.lda, param.ldc, param.n);
    //printf("%d %d %d %d\n\n", k, lda, ldc, n);

    double buffer[MAX_SIZE_DOUBLE];
    double *a_local = buffer;
    double *c_local = buffer + n;


    //printf("BEGIN calculation\n");
    //printf("%d %d %d %d\n", param.k, param.lda, param.ldc, param.n);
    //printf("%d %d %d %d\n", k, lda, ldc, n);

    for (j = 0; j < n; j++){
        //printf("j=%d\n", j);
        for (i = j; i < n; i++){
            c_local[i-j] = 0.;
        }
        for (l = 0; l < k; l++){
            get_reply = 0;
            athread_get(PE_MODE, a+(l*lda + j), a_local, \
                        sizeof(double)*(n-j), &get_reply, 0,0,0);
            while(get_reply != 1);
            temp = a_local[0];                
            if (temp > DBL_EPSILON || temp < -DBL_EPSILON){
                for (i = j; i < n; i++){
                    c_local[i - j] += temp*a_local[i - j];
                }
            } 
        }
        put_reply = 0;
        athread_put(PE_MODE, c_local, c+(j*ldc + j), \
                    sizeof(double)*(n - j), &put_reply, 0,0);
        while(put_reply != 1);
    }

    return;
}


void dsyrk_single_slave(void* ptr){
    const int my_id = athread_get_id(-1);
    if (my_id != 0) return;
    volatile unsigned int get_reply = 0, put_reply = 0;

    dsyrkParam* param = (dsyrkParam*)ptr;
    //get_reply = 0;
    //athread_get(PE_MODE, ptr, &param, 
    //            sizeof(dsyrkParam), &get_reply, 0,0,0);
    //while(get_reply != 1);

    // 在这里存在一些问题
    // athread_get param参数总是在while结束后 还没有获取完毕 
    // 当使用非-g模式时
    // 在这里借鉴了swdnn的方法，使用gload方式获取参数

    int j, l, i;
    double temp;

    int k = param->k;
    int lda = param->lda;
    int ldc = param->ldc;
    int n = param->n;
    double *a = param->a;
    double *c = param->c;

    double buffer[MAX_SIZE_DOUBLE];
    double *c_local = buffer;
    double *a_local = buffer + n;

    assert(n*k+n < MAX_SIZE_DOUBLE && "cache size");

    volatile unsigned int replygetA=0, replyputC=0;
    dma_desc dmagetA, dmaputC;
    dma_set_op(&dmagetA, DMA_GET);
    dma_set_mode(&dmagetA, PE_MODE);
    dma_set_reply(&dmagetA, &replygetA);
    dma_set_op(&dmaputC, DMA_PUT);
    dma_set_mode(&dmaputC, PE_MODE);
    dma_set_reply(&dmaputC, &replyputC);
    
    dma_set_size(&dmagetA, sizeof(double)*n*k);
    dma_set_bsize(&dmagetA, sizeof(double)*n);
    dma_set_stepsize(&dmagetA, sizeof(double)*(lda-n));    
    dma(dmagetA, (long)a, (long)a_local);
    dma_wait(&replygetA, 1);
    replygetA = 0;

    dma_set_bsize(&dmaputC, 0);
    dma_set_stepsize(&dmaputC, 0);    

    for (j = 0; j < n; j++){
        for (i = j; i < n; i++){
            c_local[i-j] = 0.;
        }
        for (l = 0; l < k; l++){
            temp = a_local[l*n+j];                
            if (temp != 0){
                for (i = j; i < n; i++){
                    c_local[i - j] += temp*a_local[l*n + i];
                }
            } 
        }
        dma_set_size(&dmaputC, sizeof(double)*(n-j));
        dma(dmaputC, (long)(c+(j*ldc + j)), (long)c_local);
        dma_wait(&replyputC, 1);
        replyputC = 0;
    }

    return;
}