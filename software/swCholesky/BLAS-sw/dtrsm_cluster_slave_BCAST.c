#include <stdlib.h>
#include <stdio.h>
#include <slave.h>
#include <assert.h>
#include <simd.h>
#include <dma.h>
#include "swBLAS_def.h"

#define DMA_A 1
#define DMA_B_IN 1
#define DMA_B_OUT 1
#define CAL 1

// TIPS : BCAST_MODE should be used with athread_syn

void dtrsm_cluster_slave(void* ptr){
    const int my_id = athread_get_id(-1);
    int k, k2, j, i;

    dtrsmParam* param = (dtrsmParam*)ptr;
    const int m = param->m;
    const int n = param->n;
    const int lda = param->lda;
    const int ldb = param->ldb;
    const double* a = param->a;
    double* b = param->b;

    int blkM = (m-1)/64+1;
    int curblkM = MAX(0, MIN(blkM, m - my_id*blkM));
    //assert(curblkM == blkM);

    volatile unsigned int replygetA=0, replygetB=0, replyputB=0;
    dma_desc dmagetA, dmagetB, dmaputB;
    dma_set_op(&dmagetA, DMA_GET);
    dma_set_mask(&dmagetA , 0xff);
    dma_set_mode(&dmagetA, BCAST_MODE);
    //dma_set_mode(&dmagetA, PE_MODE);
    dma_set_reply(&dmagetA, &replygetA);
    dma_set_bsize(&dmagetA, 0);
    dma_set_stepsize(&dmagetA, 0);

    dma_set_op(&dmagetB, DMA_GET);
    dma_set_mode(&dmagetB, PE_MODE);
    dma_set_reply(&dmagetB, &replygetB);

    dma_set_op(&dmaputB, DMA_PUT);
    dma_set_mode(&dmaputB, PE_MODE);
    dma_set_reply(&dmaputB, &replyputB);

    doublev4 vdbl, vdbl2;
    int offset, offsetL, offsetR;
    double *start;
    double *a_local; 
    double *b_local;
    int double_buffer_flag = 0;

    //doublev4 MEM[MAX_SIZE_DOUBLE/4];
    doublev4 MEM[60*1024/8/4];
    b_local = (double*)MEM;
    a_local = b_local + n*blkM;

        //a_local = (double*)((doublev4*)ldm_malloc(sizeof(double)*n*2));
#if DMA_A
        k = 0;
        start = &a[BLAS_OFFSET(k,k,lda)];
        dma_set_size(&dmagetA, sizeof(double)*(n-k));
        if (my_id == 0){
            dma(dmagetA, (long)start, (long)(a_local+k+(1-double_buffer_flag)*n) );
            dma_wait(&replygetA, 1);
        }
        athread_syn(ARRAY_SCOPE,0xffff);
        replygetA = 0;
        double_buffer_flag = 1;
#endif



#if DMA_B_IN
    //if (curblkM > 0){
        //b_local = (double*)((doublev4*)ldm_malloc(sizeof(double)*n*blkM));
        start = &b[my_id*blkM];
        dma_set_size(&dmagetB, sizeof(double)*n*curblkM);
        dma_set_bsize(&dmagetB, sizeof(double)*curblkM);
        dma_set_stepsize(&dmagetB, sizeof(double)*(ldb-curblkM));
        dma(dmagetB, (long)start, (long)b_local);
        dma_wait(&replygetB, 1);
        replygetB = 0;
    //}
#endif

    for (k = 0; k < n-1; k++){
#if DMA_A
        k2 = k+1;
        //if(k2 < n){
            start = &a[BLAS_OFFSET(k2,k2,lda)];
            dma_set_size(&dmagetA, sizeof(double)*(n-k2));
            if(my_id==0) dma(dmagetA, (long)start, (long)(a_local+k2+(1-double_buffer_flag)*n) );
        //}
#endif

#if CAL
        double temp = 1./a_local[k + double_buffer_flag*n];
        doublev4 temp_v4 = temp;
        
        for(i = 0; i < curblkM; i+=4){
            //((doublev4*)(&b_local[BLAS_OFFSET(i<<2,k,curblkM)]))[0] = \
                ((doublev4*)(&b_local[BLAS_OFFSET(i<<2,k,curblkM)]))[0] * temp_v4;
            offset = BLAS_OFFSET(i,k,curblkM);
            simd_load(vdbl, &b_local[offset]);
            vdbl *= temp_v4;
            simd_store(vdbl, &b_local[offset]);
        }

        for (j = k+1; j < n; j++){
            temp = a_local[j + double_buffer_flag*n];
            temp_v4 = temp;
            
            for(i = 0; i < curblkM; i+=4){
                //((doublev4*)(&b_local[BLAS_OFFSET(i<<2,j,curblkM)]))[0] = \
                    ((doublev4*)(&b_local[BLAS_OFFSET(i<<2,j,curblkM)]))[0] - \
                    temp_v4 * ((doublev4*)(&b_local[BLAS_OFFSET(i<<2,k,curblkM)]))[0];
                offsetL = BLAS_OFFSET(i,j,curblkM);
                offsetR = BLAS_OFFSET(i,k,curblkM);
                simd_load(vdbl,  &b_local[offsetL]);
                simd_load(vdbl2, &b_local[offsetR]);
                vdbl -= vdbl2 * temp_v4;
                simd_store(vdbl, &b_local[BLAS_OFFSET(i,j,curblkM)]);
            }
            
        }
#endif

#if DMA_A
        //if(k2 < n){
            if(my_id==0) dma_wait(&replygetA, 1);
            athread_syn(ARRAY_SCOPE,0xffff);
            replygetA = 0;
            double_buffer_flag = 1-double_buffer_flag;
        //} 
#endif
    }

#if 1
        k = n-1;
        double temp = 1./a_local[k + double_buffer_flag*n];
        doublev4 temp_v4 = temp;
        
        for(i = 0; i < curblkM; i+=4){
            //((doublev4*)(&b_local[BLAS_OFFSET(i<<2,k,curblkM)]))[0] = \
                ((doublev4*)(&b_local[BLAS_OFFSET(i<<2,k,curblkM)]))[0] * temp_v4;
            simd_load(vdbl, &b_local[BLAS_OFFSET(i,k,curblkM)]);
            vdbl *= temp_v4;
            simd_store(vdbl, &b_local[BLAS_OFFSET(i,k,curblkM)]);
        }

        for (j = k+1; j < n; j++){
            temp = a_local[j + double_buffer_flag*n];
            temp_v4 = temp;
            
            for(i = 0; i < curblkM; i+=4){
                //((doublev4*)(&b_local[BLAS_OFFSET(i<<2,j,curblkM)]))[0] = \
                    ((doublev4*)(&b_local[BLAS_OFFSET(i<<2,j,curblkM)]))[0] - \
                    temp_v4 * ((doublev4*)(&b_local[BLAS_OFFSET(i<<2,k,curblkM)]))[0];
                simd_load(vdbl, &b_local[BLAS_OFFSET(i,j,curblkM)]);
                simd_load(vdbl2, &b_local[BLAS_OFFSET(i,k,curblkM)]);
                vdbl -= vdbl2 * temp_v4;
                simd_store(vdbl, &b_local[BLAS_OFFSET(i,j,curblkM)]);
            }
            
        }
#endif

#if DMA_B_OUT
    //if (curblkM > 0){
        start = &b[my_id*blkM];
        dma_set_size(&dmaputB, sizeof(double)*n*curblkM);
        dma_set_bsize(&dmaputB, sizeof(double)*curblkM);
        dma_set_stepsize(&dmaputB, sizeof(double)*(ldb-curblkM));
        dma(dmaputB, (long)start, (long)b_local);
        dma_wait(&replyputB, 1);
        replyputB = 0;

        //ldm_free(b_local, sizeof(double)*n*blkM);
    //}
#endif

        //ldm_free(a_local, sizeof(double)*n*2);

    return;


}
