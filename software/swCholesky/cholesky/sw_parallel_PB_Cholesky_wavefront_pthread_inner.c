#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "../common/sw_Reach.h"
#include "sw_parallel_PB_Cholesky_wavefront.h"

#define __DEBUG_PRINTF__  
#ifdef __DEBUG_PRINTF__  
#define PRINT(format,...) printf("Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)  
#else  
#define PRINT(format,...)  
#endif  


void* sw_cholesky_left_par_waveFront_pthread_inner(void* ptr){
    PRINT("THREAD START inside pthread_function\n");
#if 0
    //int ret = 0;
    //ret = athread_init();
    const WavefrontInnerParam *param = (WavefrontInnerParam*)ptr;
    int n           = param->n;
    int *c          = param->c;
    int *r          = param->r;
    double *values  = param->values;

    size_t *lC      = param->lC;
    double *lValues = param->lValues;
    int *lR         = param->lR;
    size_t *Li_ptr  = param->Li_ptr;

    int *blockSet   = param->blockSet;
    int supNo       = param->supNo;

    int *prunePtr   = param->prunePtr;
    int *pruneSet   = param->pruneSet;

    //int nLevels     = param->nLevels;
    //int *levelPtr   = param->levelPtr;
    //int *levelSet   = param->levelSet;

    int super_max   = param->super_max;
    int col_max     = param->col_max;

    const int pruneBgn      = param->pruneBgn;
    const int pruneEnd      = param->pruneEnd;
    const int nSupR         = param->nSupR;
    const int curCol        = param->curCol;
    const int supWdt        = param->supWdt;
    const int* map          = param->map;
    double* cur             = param->cur;
    pthread_mutex_t* lock = param->lock;

//**********************************************************************************************************************//
 char Trana,Tranb;
 char Tran, UPLO;
 char SIDE, DIAG;

 //int *xi;
 //int *map ;
 double *contribs;
 double one = 1.0, zero = 0.0;
//**********************************************************************************************************************//
    if (pruneEnd-pruneBgn>0)
        contribs = (double*)malloc(sizeof(double)* super_max * col_max);

    int i,j;
     for (i = pruneBgn; i < pruneEnd; ++i) {
         printf("i=%d\n", i);
      int lSN = pruneSet[i];

     int nSupRs = 0;
     int cSN = blockSet[lSN];//first col of current SN
     int cNSN = blockSet[lSN+1];//first col of Next SN
     //assert(cSN < n);
     //assert(cNSN < n);
     size_t Li_ptr_cNSN = Li_ptr[cNSN];
     size_t Li_ptr_cSN = Li_ptr[cSN];
     //assert(Li_ptr_cNSN < n);
     //assert(Li_ptr_cSN < n);
     int nSNRCur=Li_ptr_cNSN-Li_ptr_cSN;
     int  supWdts=cNSN-cSN;//The width of current src SN
     int lb=0,  ub=0;
     int sw=1;
     for (j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
      //finding the overlap between curCol and curCol+supWdt in the src col
      if (lR[j] >= curCol && sw) {
       //src*transpose(row lR[j])
       lb=j-Li_ptr_cSN;
       sw=0;
      }
      if(lR[j] < curCol+supWdt && !sw){
       ub=j-Li_ptr_cSN;
      }
       if(lR[j] >= curCol + supWdt )
        break;
     }
     nSupRs=Li_ptr_cNSN-Li_ptr_cSN-lb;
     int ndrow1=ub-lb+1;
     int ndrow3 = nSupRs-ndrow1;
     double* src=&lValues[lC[cSN]+lb];//first element of src supernode starting from row lb
     double *srcL = &lValues[lC[cSN]+ub+1];
    #ifdef XMATH
      UPLO = 'L';
      Tran = 'N';
      PRINT("inner dsyrk_\n");
      PRINT("## dsyrk param :: n=%d, k=%d\n", ndrow1, supWdts);
      PRINT("## dsyrk param :: alpha=%f, beta=%f\n", one, zero);
      PRINT("## dsyrk param :: lda=%d, ldc=%d\n", nSNRCur, nSupRs);
      PRINT("## dsyrk param :: a=%p, c=%p\n", src, contribs);
     dsyrk_(&UPLO,&Tran,&ndrow1,&supWdts,&one,src,&nSNRCur,&zero,
                            contribs,&nSupRs);
      PRINT("inner dsyrk_ end\n");
    #endif
     if(ndrow3>0){
    #ifdef XMATH
       Trana = 'N'; 
       Tranb = 'C'; 
      PRINT("inner dgemm_\n");
      PRINT("## dgemm param :: m=%d, n=%d, k=%d\n", ndrow3, ndrow1, supWdts);
      PRINT("## dgemm param :: alpha=%f, beta=%f\n", one, zero);
      PRINT("## dgemm param :: lda=%d, ldb=%d, ldc=%d\n", nSNRCur, nSNRCur, nSupRs);
      PRINT("## dgemm param :: a=%p, b=%p, c=%p\n", srcL, src, contribs+ndrow1);
      dgemm_(&Trana, &Tranb,&ndrow3,&ndrow1,&supWdts,&one,srcL,&nSNRCur,
                               src,&nSNRCur,&zero,contribs+ndrow1,&nSupRs );
      PRINT("inner dgemm_ end\n");
    #endif

     }
     PRINT("## BEGIN Copy\n");
     //copying contrib to L
     pthread_mutex_lock(lock);
     int i2, j2;
     for (i2 = 0; i2 < ndrow1; ++i2) {//Copy contribs to L
      int col=map[lR[Li_ptr_cSN+i2+lb]];//col in the SN
      for (j2 = i2; j2 < nSupRs ; ++j2) {
       int cRow= lR[Li_ptr_cSN+j2+lb];//corresponding row in SN
       //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
       cur[col*nSupR+map[cRow]] -= contribs[i2*nSupRs+j2];
      }
     }
     pthread_mutex_unlock(lock);
     PRINT("## FINISH Copy\n");
    }//end top

    if (pruneEnd-pruneBgn>0)
        free(contribs);
#endif
    PRINT("THREAD FINISH inside pthread_function\n");

    return NULL;
}


