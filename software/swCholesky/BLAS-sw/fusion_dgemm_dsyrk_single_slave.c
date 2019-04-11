#include <stdlib.h>
#include <stdio.h>
#include <float.h>
//#include <assert.h>
#include <slave.h>
#include <dma.h>
#include "swBLAS_def.h"


void fusion_dgemm_dsyrk_single_slave(void* ptr){
    const int my_id = athread_get_id(-1);
#ifdef DEBUG_SINGLE
    if (my_id != 0) return;
    fusionDgemmDsyrkParam* param = (fusionDgemmDsyrkParam*)ptr;
#else
    fusionTPParam* TPParam = (fusionTPParam*)ptr;
    const int TPnum = TPParam->num;
    if (my_id >= TPnum) return;
#endif

    fusionDgemmDsyrkParam *param = &(TPParam->Q)[my_id];
    //printf("%d\n", my_id);
    //printf("size = %d\n", sizeof(fusionDgemmDsyrkParam) );
    //printf("Q = %p\n", (TPParam->Q)+my_id);
/*
    fusionDgemmDsyrkParam param;
    volatile unsigned int replygetP=0;
    dma_desc dmagetP;
    dma_set_op(&dmagetP, DMA_GET);
    dma_set_mode(&dmagetP, PE_MODE);
    dma_set_reply(&dmagetP, &replygetP);
    dma_set_size(&dmagetP, sizeof(fusionDgemmDsyrkParam));
    dma_set_bsize(&dmagetP, 0);
    dma_set_stepsize(&dmagetP, 0);
    dma(dmagetP, (long)((TPParam->Q)+my_id), (long)(&param)); 
    dma_wait(&replygetP, 1);
*/
    int lda = param->lda;
    int ldc = param->ldc;
    int ndrow1 = param->ndrow1;
    int ndrow3 = param->ndrow3;
    int wdt = param->wdt;
    double* a = param->a;
    double* c = param->c;
    int lda_local = ndrow1+ndrow3;

    int j, l, i;
    double temp;

    double buffer[MAX_SIZE_DOUBLE];
    double* c_local = buffer;
    double* a_local = c_local + ndrow1+ndrow3;

    //printf("my_id is %d :: ndrow1=%d, ndrow3=%d, wdt=%d\n", my_id,ndrow1,ndrow3,wdt);
    //printf("a = %p, c = %p\n", a,c);
    //assert((ndrow1+ndrow3)*(wdt+1) <= MAX_SIZE_DOUBLE && "cache size");

    volatile unsigned int replygetA=0, replyputC=0;
    dma_desc dmagetA, dmaputC;
    dma_set_op(&dmagetA, DMA_GET);
    dma_set_mode(&dmagetA, PE_MODE);
    dma_set_reply(&dmagetA, &replygetA);
    dma_set_op(&dmaputC, DMA_PUT);
    dma_set_mode(&dmaputC, PE_MODE);
    dma_set_reply(&dmaputC, &replyputC);

    dma_set_size(&dmagetA, sizeof(double)*lda_local*wdt);
    dma_set_bsize(&dmagetA, sizeof(double)*lda_local);
    dma_set_stepsize(&dmagetA, sizeof(double)*(lda-lda_local));  

    dma(dmagetA, (long)a, (long)a_local);
    dma_wait(&replygetA, 1);
    replygetA = 0;

    dma_set_bsize(&dmaputC, 0);
    dma_set_stepsize(&dmaputC, 0);    


    for (j = 0; j < ndrow1; j++){
        for (i = j; i < lda_local; i++){
            temp = 0.;
            for (l = 0; l < wdt; l++){
                temp += a_local[l*lda_local + i] * a_local[l*lda_local + j];
            }
            c_local[i] = temp;
        }
        dma_set_size(&dmaputC, sizeof(double)*(lda_local-j));
        dma(dmaputC, (long)(c+j*ldc+j), (long)(c_local+j));
        dma_wait(&replyputC, 1);
        replyputC = 0;
    }
    
    return;
}