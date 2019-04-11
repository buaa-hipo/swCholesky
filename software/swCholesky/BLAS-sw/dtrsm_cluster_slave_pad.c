#include <stdlib.h>
#include <stdio.h>
#include <slave.h>
#include <simd.h>
#include <dma.h>
#include "swBLAS_def.h"

void dtrsm_pad_b(void* ptr){
    const int my_id = athread_get_id(-1);
    dtrsmPadParam* param = (dtrsmPadParam*)ptr;
    const int m = param->m;
    const int n = param->n;
    const int ldb = param->ldb;
    const int blkM = param->blkM;
    double* b = param->b;
    double* bp = param->bp;
    const Ms = m/blkM*blkM;
    const Me = (m+blkM-1)/blkM*blkM;
    
    if (Ms >= Me) return;

    if ( blkM%256 != 0 ){
        printf("dtrsm_pad_b error : blkM=%d\n", blkM);
        return;
    }

    //if (my_id == 0) printf("m=%d, n=%d, blkM=%d, Ms=%d, Me=%d\n", m,n,blkM,Ms,Me);
    //if (my_id != 0) return;

    volatile int reply_src = 0, reply_dst = 0;
    dma_desc dma_src, dma_dst;
    dma_set_op(&dma_src, DMA_GET);
    dma_set_mode(&dma_src, PE_MODE);
    dma_set_reply(&dma_src, &reply_src);
    dma_set_size(&dma_src, sizeof(double)*(m-Ms));
    dma_set_bsize(&dma_src, 0);
    dma_set_stepsize(&dma_src, 0);

    dma_set_op(&dma_dst, DMA_PUT);
    dma_set_mode(&dma_dst, PE_MODE);
    dma_set_reply(&dma_dst, &reply_dst);
    dma_set_size(&dma_dst, sizeof(double)*(Me-Ms));
    dma_set_bsize(&dma_dst, 0);
    dma_set_stepsize(&dma_dst, 0);

    double* buf = (double*)ldm_malloc(sizeof(double)*blkM);
    double* start;

    int i, j;
    for (j = my_id; j < n; j+=64){
        //load
        start = &b[j*ldb + Ms];
        dma(dma_src, (long)start, (long)buf);
        dma_wait(&reply_src, 1);
        reply_src = 0;
        for (i = m-Ms; i < Me-Ms; i++){
            buf[i] = 0.0;
        }

        //store
        start = &bp[j*blkM];
        dma(dma_dst, (long)start, (long)buf);
        dma_wait(&reply_dst, 1);
        reply_dst = 0;
    }

    ldm_free(buf, sizeof(double)*blkM);

    return;
}

void dtrsm_depad_b(void* ptr){
    const int my_id = athread_get_id(-1);
    dtrsmPadParam* param = (dtrsmPadParam*)ptr;
    const int m = param->m;
    const int n = param->n;
    const int ldb = param->ldb;
    const int blkM = param->blkM;
    double* b = param->b;
    double* bp = param->bp;
    const Ms = m/blkM*blkM;
    const Me = (m+blkM-1)/blkM*blkM;

    if (Ms >= Me) return;

    if ( blkM%256 != 0 ){
        printf("dtrsm_pad_b error : blkM=%d\n", blkM);
        return;
    }

    //if (my_id == 0) printf("m=%d, n=%d, blkM=%d, Ms=%d, Me=%d\n", m,n,blkM,Ms,Me);
    //if (my_id != 0) return;

    volatile int reply_src = 0, reply_dst = 0;
    dma_desc dma_src, dma_dst;
    dma_set_op(&dma_src, DMA_GET);
    dma_set_mode(&dma_src, PE_MODE);
    dma_set_reply(&dma_src, &reply_src);
    dma_set_size(&dma_src, sizeof(double)*(Me-Ms));
    dma_set_bsize(&dma_src, 0);
    dma_set_stepsize(&dma_src, 0);

    dma_set_op(&dma_dst, DMA_PUT);
    dma_set_mode(&dma_dst, PE_MODE);
    dma_set_reply(&dma_dst, &reply_dst);
    dma_set_size(&dma_dst, sizeof(double)*(m-Ms));
    dma_set_bsize(&dma_dst, 0);
    dma_set_stepsize(&dma_dst, 0);

    double* buf = (double*)ldm_malloc(sizeof(double)*blkM);
    double* start;

    int i, j;
    for (j = my_id; j < n; j+=64){
        //load
        start = &bp[j*blkM];
        dma(dma_src, (long)start, (long)buf);
        dma_wait(&reply_src, 1);
        reply_src = 0;

        //store
        start = &b[j*ldb+Ms];
        dma(dma_dst, (long)start, (long)buf);
        dma_wait(&reply_dst, 1);
        reply_dst = 0;
    }

    ldm_free(buf, sizeof(double)*blkM);

    return;
}
