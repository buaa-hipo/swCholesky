
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <athread.h>
#include <assert.h>
#include <sys/time.h>
#include "../common/sw_Reach.h"
#include "sw_parallel_PB_Cholesky_05.h"
#include "../BLAS-sw/swBLAS_def.h"
#include "../BLAS-sw/swBLAS_fun.h"
#include "tripleQueue.h"

#define TIMING
#undef TIMING
#define TIMING1
#undef TIMING1
#define BLASTIMING
#undef BLASTIMING


#define MKL_INT int
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

//#define __DEBUG_PRINTF__  
#ifdef __DEBUG_PRINTF__  
#define PRINT(format,...) printf("Line: %05d: "format"\n", __LINE__, ##__VA_ARGS__)  
#else  
#define PRINT(format,...)  
#endif  

extern int isFiniteNumber(double d);


void sw_cholesky_left_par_05_pthread_inner(void *ptr){
    PRINT("THREAD START inside pthread_function\n");
    int ret = 0;
    ret = athread_init();

#ifdef TRIPLEQUEUE
    //fusionQueue Q[64];
    //fusionDgemmDsyrkParam S[64];
    fusionQueue *Q = (fusionQueue*)malloc(sizeof(fusionQueue)*64);
    fusionDgemmDsyrkParam *S = (fusionDgemmDsyrkParam*)malloc(sizeof(fusionDgemmDsyrkParam)*64);
    int cntQ = 0;
#endif

#if 1
    const InnerParam *param = (InnerParam*)ptr;
    int n           = param->n;
    int *c          = param->c;
    int *r          = param->r;
    double *values  = param->values;

    size_t *lC      = param->lC;
    double *lValues = param->lValues;
    int *Ls         = param->Ls;
    size_t *Li_ptr  = param->Li_ptr;

    int *blockSet   = param->blockSet;
    int supNo       = param->supNo;

  #ifndef PRUNE
    int *aTree      = param->aTree;
    int *cT         = param->cT;
    int *rT         = param->rT;
    int *col2Sup    = param->col2Sup;
  #else
    int *prunePtr   = param->prunePtr;
    int *pruneSet   = param->pruneSet;
  #endif

    //int nLevels     = param->nLevels;
    //int *levelPtr   = param->levelPtr;
    //int *levelSet   = param->levelSet;

    int nPar        = param->nPar;
    int *parPtr     = param->parPtr;
    int *partition  = param->partition;

    int super_max   = param->super_max;
    int col_max     = param->col_max;
    double *nodCost = param->nodCost;

    const int levelBgn      = param->levelBgn;
    const int levelEnd      = param->levelEnd;

//**********************************************************************************************************************//

 //
 // For timing using BLAS
 //
 int top = 0;
 // int *xi = new int[2*supNo]();
 //int col_max = n;
 int *map, *xi;
 double *contribs;
 MKL_INT info = 0;
 int thth=0;
 double one, zero;
 one = 1.0;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
 zero = 0.;     /* BETA for *syrk, *herk, and *gemm */
 struct timeval start,end, startin,endin;
 double elapsed_seconds;
#ifdef BLASTIMING
 struct timeval startBlas, endBlas;
 double blasTime = 0;
 double syrkTime = 0, gemmTime = 0, potrfTime = 0, trsmTime = 0;
#endif
 double duration4 = 0 ,duration3 = 0, duration2=0, duration1=0;
 char Trana,Tranb;
 char Tran, UPLO;
 char SIDE, DIAG;

   int j1, k1, i1;
   int i, j, cnt;

   for (j1 = levelBgn; j1 < levelEnd; ++j1) {
    map = (int*)malloc(sizeof(int)* n);
    contribs = (double*)malloc(sizeof(double)* super_max * col_max);
    xi = (int*)malloc(sizeof(int)* 2*supNo);

#ifdef TIMING1
    gettimeofday(&startin, NULL);
#endif

    PRINT("\t\t ITERATION :: j1 = %d #k1 From:%d To:%d\n", j1, parPtr[j1], parPtr[j1+1]);
    for (k1 = parPtr[j1]; k1 < parPtr[j1 + 1]; ++k1) {
     int s = partition[k1] + 1;
     int curCol = s != 0 ? blockSet[s - 1] : 0;
     int nxtCol = blockSet[s];
     MKL_INT supWdt = nxtCol - curCol;
     MKL_INT nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
     for (i = Li_ptr[curCol], cnt = 0; i < Li_ptr[nxtCol]; ++i) {
      map[Ls[i]] = cnt++;//mapping L rows position to actual row idx
     }
     //copy the columns from A to L
     ///only copy lower triangular part
     for (i = curCol; i < nxtCol; ++i) {//Copy A to L
      int pad = i - curCol;
      for (j = c[i]; j < c[i + 1]; ++j) {
       // if(r[j]>=i)//does not need to save upper part.
       lValues[lC[i] + map[r[j]]] = values[j];
       //   else
       //      printf("dddd\n");
      }
     }
     double *src, *cur = &lValues[lC[curCol]];//pointing to first element of the current supernode

#ifndef PRUNE
     #ifdef TRIPLEQUEUE
     cntQ=0;
     #endif
     top = sw_ereach_sn(supNo,cT,rT,curCol,nxtCol,col2Sup, aTree,xi,xi+supNo);
     assert(top>=0);
     PRINT("\t\t\t ITERATION :: k1 = %d #i From:%d To:%d\n", k1, top, supNo);
     for(i = top; i < supNo; ++i){
      int lSN = xi[i];

  #if DEBUG
      if(xi[top++] != lSN)
                         printf("fail");
  #endif

     int isLast = 0;
     if (i == supNo - 1) isLast = 1;
#else
     #ifdef TRIPLEQUEUE
     cntQ=0;
     #endif
     PRINT("\t\t\t ITERATION :: k1 = %d #i PRUNE From:%d To:%d\n", k1, prunePtr[s-1], prunePtr[s]);
     for (i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
      int lSN = pruneSet[i];

     int isLast = 0;
     if (i == prunePtr[s] - 1) isLast = 1;
#endif
      MKL_INT nSupRs = 0;
      int cSN = blockSet[lSN];//first col of current SN
      int cNSN = blockSet[lSN + 1];//first col of Next SN
      size_t Li_ptr_cNSN = Li_ptr[cNSN];
      size_t Li_ptr_cSN = Li_ptr[cSN];
      MKL_INT nSNRCur = Li_ptr_cNSN - Li_ptr_cSN;
      MKL_INT supWdts = cNSN - cSN;//The width of current src SN
      int lb = 0, ub = 0;
      //bool sw = true;
      int sw = 1;
      for (j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
       //finding the overlap between curCol and curCol+supWdt in the src col
       if (Ls[j] >= curCol && sw) { ///???????why Ls?
        //src*transpose(row Ls[j])
        lb = j - Li_ptr_cSN;
        sw = 0;
       }
       if (Ls[j] < curCol + supWdt && !sw) {
        ub = j - Li_ptr_cSN;
       }
       if(Ls[j] >= curCol + supWdt )
        break;
      }
      nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
      MKL_INT ndrow1 = ub - lb + 1;
      MKL_INT ndrow3 = nSupRs - ndrow1;
      src = &lValues[lC[cSN] +
                     lb];//first element of src supernode starting from row lb
      double *srcL = &lValues[lC[cSN] + ub + 1];
#ifdef BLASTIMING
      gettimeofday(&startBlas, NULL);
#endif
      PRINT("inner step : 1, %d\n", i);

      //int mz;
      //for (mz = 0; mz < super_max*col_max; mz++)
      //for (mz = 0; mz < ndrow1*nSupRs; mz++)
      //  contribs[mz] = 0.0;

#ifndef TRIPLEQUEUE

    #ifdef XMATH
      //printf("inner dsyrk_\n");
      UPLO = 'L';
      Tran = 'N';
      dsyrk_(&UPLO,&Tran,&ndrow1,&supWdts,&one,src,&nSNRCur,&zero,
                             contribs,&nSupRs);
      //printf("inner dsyrk_ end\n");
    #endif

    #ifdef MYBLAS
      PRINT("inner dsyrk_\n");
      PRINT("## dsyrk param :: n=%d, k=%d\n", ndrow1, supWdts);
      PRINT("## dsyrk param :: alpha=%f, beta=%f\n", one, zero);
      PRINT("## dsyrk param :: lda=%d, ldc=%d\n", nSNRCur, nSupRs);
      PRINT("## dsyrk param :: a=%p, c=%p\n", src, contribs);
      dsyrk_cluster_master(ndrow1,supWdts,one,src,nSNRCur,zero,
                             contribs,nSupRs);
      PRINT("inner dsyrk_ end\n");
    #endif
      PRINT("inner step : 2, %d\n", i);
        #ifdef BLASTIMING
       gettimeofday(&endBlas, NULL);
       elapsed_seconds = TIME(startBlas, endBlas);
       syrkTime += elapsed_seconds;
       blasTime += elapsed_seconds;
        #endif

      if (ndrow3 > 0) {
        #ifdef BLASTIMING
      gettimeofday(&startBlas, NULL);
        #endif
    #ifdef XMATH
      //printf("inner dgemm_\n");
       Trana = 'N'; 
       Tranb = 'C'; 
       dgemm_(&Trana,&Tranb,&ndrow3,&ndrow1,&supWdts,&one,srcL,&nSNRCur,
                                src,&nSNRCur,&zero,contribs+ndrow1,&nSupRs );
      //printf("inner dgemm_ end\n");
    #endif

    #ifdef MYBLAS
      PRINT("inner dgemm_\n");
      PRINT("## dgemm param :: m=%d, n=%d, k=%d\n", ndrow3, ndrow1, supWdts);
      PRINT("## dgemm param :: alpha=%f, beta=%f\n", one, zero);
      PRINT("## dgemm param :: lda=%d, ldb=%d, ldc=%d\n", nSNRCur, nSNRCur, nSupRs);
      PRINT("## dgemm param :: a=%p, b=%p, c=%p\n", srcL, src, contribs+ndrow1);
       dgemm_cluster_master(ndrow3,ndrow1,supWdts,one,srcL,nSNRCur,
                                src,nSNRCur,zero,contribs+ndrow1,nSupRs );
      PRINT("inner dgemm_ end\n");
    #endif
        #ifdef BLASTIMING
       gettimeofday(&endBlas, NULL);
       elapsed_seconds = TIME(startBlas, endBlas);
       gemmTime += elapsed_seconds;
       blasTime += elapsed_seconds;
        #endif
      }// end of ndrow3
      
      PRINT("inner step : 3, %d\n", i);

      //copying contrib to L
      int i2, j2;
      for (i2 = 0; i2 < ndrow1; ++i2) {//Copy contribs to L
       int col = map[Ls[Li_ptr_cSN + i2 + lb]];//col in the SN
       for (j2 = i2; j2 < nSupRs; ++j2) {
        int cRow = Ls[Li_ptr_cSN + j2 + lb];//corresponding row in SN
        //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
      #ifdef DEBUG_CHECK
        assert(!isnan(cur[col * nSupR + map[cRow]]));
        assert(!isnan(contribs[i2 * nSupRs + j2]));
        assert(isFiniteNumber(cur[col * nSupR + map[cRow]]));
        assert(isFiniteNumber(contribs[i2 * nSupRs + j2]));
      #endif
        cur[col * nSupR + map[cRow]] -= contribs[i2 * nSupRs + j2];
       }
      }
      PRINT("inner step : 4, %d\n", i);
#else
    PRINT("inner queue : 0, %d\n", i);
    int retval;
    retval = addTask(Q, &cntQ, ndrow1, ndrow3, nSupRs, supWdts, nSNRCur, \
                    Li_ptr_cSN, lb, ub, src, srcL, contribs);
    PRINT("retval = %d, cntQ = %d, isLast = %d\n", retval, cntQ, isLast);
    if ((retval == -1) || isLast) {
        PRINT("EXECUTE QUEUE, size=%d, isLast=%d#################\n", cntQ, isLast);
        execQueue(Q, S, &cntQ, nSupR, map, Ls, cur);
    }
    if (retval == 0 || retval == 2) {
      //copying contrib to L
      int i2, j2;
      for (i2 = 0; i2 < ndrow1; ++i2) {//Copy contribs to L
       int col = map[Ls[Li_ptr_cSN + i2 + lb]];//col in the SN
       for (j2 = i2; j2 < nSupRs; ++j2) {
        int cRow = Ls[Li_ptr_cSN + j2 + lb];//corresponding row in SN
        //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
      #ifdef DEBUG_CHECK
        assert(!isnan(cur[col * nSupR + map[cRow]]));
        assert(!isnan(contribs[i2 * nSupRs + j2]));
        assert(isFiniteNumber(cur[col * nSupR + map[cRow]]));
        assert(isFiniteNumber(contribs[i2 * nSupRs + j2]));
      #endif
        cur[col * nSupR + map[cRow]] -= contribs[i2 * nSupRs + j2];
       }
      }
    }
    PRINT("inner queue : 1, %d\n", i);
#endif



#ifndef PRUNE
     } /// end of pruneSet
#else
     } /// end of none pruneSet
#endif

#ifdef BLASTIMING
     gettimeofday(&startBlas, NULL);
#endif
#ifdef XMATH
     //printf("inner dpotrf \n");
     char UPLO = 'L';
     dpotrf_(&UPLO,&supWdt,cur,&nSupR,&info);
     //printf("inner dpotrf end\n");
     if(info!=0){
         printf("dpotrf :: nonpositive definite matrix, info=%d\n", info);
         break;
     }
#endif
#ifdef MYBLAS
     PRINT("inner dpotrf \n");
     //Cholesky_col(nSupR,supWdt,cur);
     PRINT("## dpotrf param :: n=%d, lda=%d\n", supWdt, nSupR);
     PRINT("## dpotrf param :: a=%p\n", cur);
     dpotrf_cluster_master(supWdt,cur,nSupR,&info);
     if(info!=0){
         printf("dpotrf :: nonpositive definite matrix, info=%d\n", info);
         break;
         exit(0);
     }
     PRINT("inner dpotrf end\n");
#endif
#ifdef BLASTIMING
       gettimeofday(&endBlas, NULL);
       elapsed_seconds = TIME(startBlas, endBlas);
       potrfTime += elapsed_seconds;
       blasTime += elapsed_seconds;
#endif

#ifdef BLASTIMING
     gettimeofday(&startBlas, NULL);
#endif
     MKL_INT rowNo = nSupR - supWdt;
#ifdef XMATH
     //printf("inner dtrsm \n");
     SIDE = 'R';
     UPLO = 'L';
     Trana = 'C';
     DIAG = 'N';
     dtrsm_(&SIDE, &UPLO, &Trana, &DIAG, &rowNo, &supWdt,&one,
                        cur,&nSupR,&cur[supWdt],&nSupR);
     //printf("inner dtrsm end\n");
#endif
#ifdef MYBLAS
     PRINT("inner dtrsm \n");
     PRINT("## dtrsm param :: m=%d, n=%d, alpha=%f, lda=%d, ldb=%d\n", rowNo, supWdt, one, nSupR, nSupR);
     PRINT("## dtrsm param :: a=%p, b=%p\n", cur, &cur[supWdt]);
     dtrsm_cluster_master(rowNo, supWdt,one,
                        cur,nSupR,&cur[supWdt],nSupR);
     PRINT("inner dtrsm end\n");
#endif
#ifdef BLASTIMING
     gettimeofday(&endBlas, NULL);
     elapsed_seconds = TIME(startBlas, endBlas);
     trsmTime += elapsed_seconds;
     blasTime += elapsed_seconds;
#endif

    } //end of k1


#ifdef TIMING1
    //endin = std::chrono::system_clock::now();
    //elapsed_seconds = endin-startin;
    //duration1=elapsed_seconds.count();
    //int thth2=omp_get_thread_num();
    //std::cout<<"**"<<thth2<<" : "<<j1<<" "<<duration1<<"\n";
    gettimeofday(&endin, NULL);
    elapsed_seconds = TIME(startin, endin);
    duration1 = elapsed_seconds;
    //std::cout<<"**"<<" : "<<j1<<" "<<duration1<<"\n" << std::endl;
    printf("** : %d --> %f \n", j1, duration1);
#endif


    free(map);
    free(contribs);
    free(xi);
   }//end of j1
#endif


#ifdef TRIPLEQUEUE
    free(Q);
    free(S);
#endif
    PRINT("THREAD FINISH inside pthread_function\n");

}
