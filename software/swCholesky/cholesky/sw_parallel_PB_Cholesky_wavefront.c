#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "../common/sw_Reach.h"
#include "sw_parallel_PB_Cholesky_wavefront.h"

//#define __DEBUG_PRINTF__  
#ifdef __DEBUG_PRINTF__  
#define PRINT(format,...) printf("Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)  
#else  
#define PRINT(format,...)  
#endif  

//#undef PTHREAD

#define NUM_CG 1
extern void* sw_cholesky_left_par_waveFront_pthread_inner(void* ptr);


int sw_cholesky_left_par_waveFront(int n, int* c, int* r, double* values,
                                 size_t *lC, int* lR, size_t* Li_ptr, double* lValues,
                          int *blockSet, int supNo, double *timing,
                          int *prunePtr, int *pruneSet,
                          int nLevels, int *levelPtr, int *levelSet,
                          int chunk, int threads,int super_max, int col_max) {
printf("####wavefront####\n");                        
printf("n = %d, supNo = %d\n", n, supNo);
printf("nLevels = %d\n", nLevels);
printf("chunk = %d, threads = %d\n", chunk, threads);
printf("super_max=%d, col_max=%d\n", super_max, col_max);


 char Trana,Tranb;
 char Tran, UPLO;
 char SIDE, DIAG;

 int top=0;
 int *xi;
 int *map ;
 double *contribs ;
 int info;
 const double one = 1.0, zero = 0.0;

#ifdef PTHREAD
    int pcnt=0;
    pthread_mutex_t lock;
    pthread_t* pthread_handler;
    pthread_handler = (pthread_t*)malloc(sizeof(pthread_t)*4);

    WavefrontInnerParam *params;
    params = (WavefrontInnerParam*)malloc(sizeof(WavefrontInnerParam)*4);
    int iii;
    for(iii = 0; iii < NUM_CG; iii++){
	    (params+iii)->n = n;
	    (params+iii)->c = c;
	    (params+iii)->r = r;
	    (params+iii)->values = values;
	    (params+iii)->lC = lC;
	    (params+iii)->lR = lR;
	    (params+iii)->Li_ptr = Li_ptr;
	    (params+iii)->lValues = lValues;
	    (params+iii)->blockSet = blockSet;
	    (params+iii)->supNo = supNo;

	    (params+iii)->prunePtr = prunePtr;
	    (params+iii)->pruneSet = pruneSet;

	    //params+iii->nLevels = nLevels;
	    //params+iii->levelPtr = levelPtr;
	    //params+iii->levelSet = levelSet;
	    (params+iii)->super_max = super_max;
	    (params+iii)->col_max = col_max;

        (params+iii)->lock = &lock;
    }
#endif


int lev, lIter, i, j, cnt;
 printf("@@START ITERATION :: lev, nLevels=%d\n", nLevels);
 for (lev = 0; lev < nLevels; ++lev) {
  
   map = (int*)malloc(sizeof(int)*n);
   contribs = (double*)malloc(sizeof(double)* super_max * col_max);
   xi = (int*)malloc(sizeof(int)* (2*n+1));

    printf("\t ITERATION :: lIter #lev From:%d To:%d\n", levelPtr[lev], levelPtr[lev+1]);
   for (lIter = levelPtr[lev]; lIter < levelPtr[lev+1]; ++lIter) {
       PRINT("lIter=%d\n", lIter);
    int s = levelSet[lIter]+1;
    int curCol = s!=0 ? blockSet[s - 1] : 0;
    int nxtCol = blockSet[s];
    int supWdt = nxtCol-curCol;
    int nSupR = Li_ptr[nxtCol]-Li_ptr[curCol];//row size of supernode
    for (i = Li_ptr[curCol],cnt=0; i < Li_ptr[nxtCol]; ++i) {
     map[lR[i]] = cnt++;//mapping L rows position to actual row idx
    }
    //copy the columns from A to L
    for (i = curCol; i < nxtCol; ++i) {//Copy A to L
     int pad=i-curCol;
     for (j = c[i]; j < c[i+1] ; ++j) {
      // if(r[j]>=i)//does not need to save upper part.
      lValues[lC[i]+map[r[j]]] = values[j];
      //   else
      //      printf("dddd\n");
     }
    }

    double *src, *cur=&lValues[lC[curCol]];//pointing to first element of the current supernode
#ifndef PTHREAD
     for (i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
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
     src=&lValues[lC[cSN]+lb];//first element of src supernode starting from row lb
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
     int i2, j2;
     for (i2 = 0; i2 < ndrow1; ++i2) {//Copy contribs to L
      int col=map[lR[Li_ptr_cSN+i2+lb]];//col in the SN
      for (j2 = i2; j2 < nSupRs ; ++j2) {
       int cRow= lR[Li_ptr_cSN+j2+lb];//corresponding row in SN
       //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
       cur[col*nSupR+map[cRow]] -= contribs[i2*nSupRs+j2];
      }
     }
     PRINT("## FINISH Copy\n");
    }//end top
#else
    printf("prepare pthread param\n");
    pthread_mutex_init(&lock, NULL);

    int pruneBgn = prunePtr[s-1];
    int pruneEnd = prunePtr[s];
    int pruneNUM = pruneEnd - pruneBgn;
    int cntcnt = (pruneNUM-1)/NUM_CG + 1;

    printf("pruneBgn = %d, pruneEnd = %d, pruneNUM = %d, cntcnt = %d\n", pruneBgn, pruneEnd, pruneNUM, cntcnt);

    for (iii = 0; iii < NUM_CG; iii++){
        (params+iii)->pruneBgn = pruneBgn + iii * cntcnt;
        (params+iii)->pruneEnd = (iii+1)*cntcnt > pruneNUM ? pruneEnd : pruneBgn+(iii+1)*cntcnt;

        (params+iii)->nSupR = nSupR;
        (params+iii)->curCol = curCol;
        (params+iii)->supWdt = supWdt;
        (params+iii)->cur = cur;
        (params+iii)->map = map;

        printf("param[%d], bgn=%d, end=%d\n", iii, (params+iii)->pruneBgn, (params+iii)->pruneEnd);
    }

    printf("pthread param finish\n");
    printf("pthread NUM_CG=%d\n", NUM_CG);

    for(iii = 0; iii < NUM_CG; iii++){
        int rc = pthread_create(&pthread_handler[iii], NULL, sw_cholesky_left_par_waveFront_pthread_inner, NULL );
        if (rc != 0)
            printf("ERROR ##FINISH pthread create%d, rc=%d\n", iii, rc);
    }
    for(iii = 0; iii < NUM_CG; iii++){
        int rc = pthread_join(pthread_handler[iii], NULL);
        if (rc != 0)
            printf("ERROR ##FINISH pthread join%d, rc=%d\n", iii, rc);
    }

    pthread_mutex_destroy(&lock);
    printf("##FINISH pthread s-1 = %d, pcnt=%d\n", s-1, pcnt);
    pcnt++;
#endif

    #ifdef XMATH
    UPLO = 'L';
     PRINT("inner dpotrf \n");
     PRINT("## dpotrf param :: n=%d, lda=%d\n", supWdt, nSupR);
     PRINT("## dpotrf param :: a=%p\n", cur);
    dpotrf_(&UPLO,&supWdt,cur,&nSupR,&info);
     if(info!=0){
         printf("dpotrf :: nonpositive definite matrix, info=%d\n", info);
         exit(0);
         break;
     }
     PRINT("inner dpotrf end\n");
    #endif

    int rowNo=nSupR-supWdt;
    #ifdef XMATH
     SIDE = 'R';
     UPLO = 'L';
     Trana = 'C';
     DIAG = 'N';
     PRINT("inner dtrsm \n");
     PRINT("## dtrsm param :: m=%d, n=%d\n", rowNo, supWdt);
     PRINT("## dtrsm param :: alpha=%f\n", one);
     PRINT("## dtrsm param :: lda=%d, ldb=%d\n", nSupR, nSupR);
     PRINT("## dtrsm param :: m=%d, n=%d, alpha=%f, lda=%d, ldb=%d\n", rowNo, supWdt, one, nSupR, nSupR);
     PRINT("## dtrsm param :: a=%p, b=%p\n", cur, &cur[supWdt]);
    dtrsm_(&SIDE, &UPLO, &Trana, &DIAG, &rowNo, &supWdt,&one,
                       cur,&nSupR,&cur[supWdt],&nSupR);
     PRINT("inner dtrsm end\n");
    #endif


   }//end lIter
   free(xi);
   free(contribs);
   free(map);

 }//end Level

#ifdef PTHREAD  
    free(params);
    free(pthread_handler);
#endif

 return 1;
}

