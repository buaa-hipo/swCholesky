
#include <stdio.h>
#include <stdlib.h>
#include <athread.h>
#include <assert.h>
#include <float.h>
#include <sys/time.h>
#include "../common/sw_Reach.h"
#include "sw_parallel_PB_Cholesky_05.h"
#include "../BLAS-sw/swBLAS_fun.h"

#define TIMING
#undef TIMING
#define TIMING1
#undef TIMING1
#define BLASTIMING
#undef BLASTIMING

#ifdef MYBLAS
    #undef MYBLAS
    #define XMATH
#endif

#ifdef PTHREAD
    #undef PTHREAD
#endif

#define MKL_INT int
#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

//#define __DEBUG_PRINTF__  
#ifdef __DEBUG_PRINTF__  
#define PRINT(format,...) printf("File: "__FILE__", Line: %05d: "format"/n", __LINE__, ##__VA_ARGS__)  
#else  
#define PRINT(format,...)  
#endif  

#define NUM_CG 4

static int isFiniteNumber(double d) {return (d<=DBL_MAX&&d>=-DBL_MAX);}
//extern void sw_cholesky_left_par_05_pthread_inner(void *ptr);

int sw_cholesky_left_par_05_check(int n, int* c, int* r, double* values,
                          size_t *lC, int* Ls, size_t* Li_ptr, double* lValues,
                          int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                          int *aTree, int *cT, int *rT, int *col2Sup,
#else
                          int *prunePtr, int *pruneSet,
#endif

                          int nLevels, int *levelPtr, int *levelSet,
                          int nPar, int *parPtr, int *partition,
                          int chunk, int threads, int super_max,
                          ///int col_max, double *nodCost=NULL
                          int col_max, double *nodCost
                          ) {
#ifndef PTHREAD                              
  athread_init();
#endif
 /*
  * For timing using BLAS
  */
 int top = 0;
 int *map, *xi;
 double *contribs;
 MKL_INT info;
 int thth=0;
 double one, zero;
 one = 1.0;    /* ALPHA for *syrk, *herk, *gemm, and *trsm */
 zero = 0.;     /* BETA for *syrk, *herk, and *gemm */
 struct timeval start,end, startin,endin;
 double elapsed_seconds;
 double duration4 = 0 ,duration3 = 0, duration2 = 0, duration1 = 0;
#ifdef BLASTIMING
 struct timeval startBlas, endBlas;
 double blasTime = 0;
 double syrkTime = 0, gemmTime = 0, potrfTime = 0, trsmTime = 0;
#endif
 
 char Trana,Tranb;
 char Tran, UPLO;
 char SIDE, DIAG;

#ifdef PTHREAD
    pthread_t pthread_handler[NUM_CG];

    int iii;
    InnerParam *params[NUM_CG];
    for(iii = 0; iii < NUM_CG; iii++){
        params[iii] = malloc(sizeof(InnerParam));
    }
    for(iii = 0; iii < NUM_CG; iii++){
	    params[iii]->n = n;
	    params[iii]->c = c;
	    params[iii]->r = r;
	    params[iii]->values = values;
	    params[iii]->lC = lC;
	    params[iii]->lValues = lValues;
	    params[iii]->Ls = Ls;
	    params[iii]->Li_ptr = Li_ptr;
	    params[iii]->blockSet = blockSet;
	    params[iii]->supNo = supNo;
	  #ifndef PRUNE
	    params[iii]->aTree = aTree;
	    params[iii]->cT = cT;
	    params[iii]->rT = rT;
	    params[iii]->col2sup = col2sup;
	  #else
	    params[iii]->prunePtr = prunePtr;
	    params[iii]->pruneSet = pruneSet;
	  #endif
	    params[iii]->nLevels = nLevels;
	    params[iii]->levelPtr = levelPtr;
	    params[iii]->levelSet = levelSet;
	    params[iii]->nPar = nPar;
	    params[iii]->parPtr = parPtr;
	    params[iii]->partition = partition;
	    params[iii]->super_max = super_max;
	    params[iii]->col_max = col_max;
	    params[iii]->nodCost = nodCost;
    }
#endif

#ifdef TIMING
   // record the whole running time, part 0
   gettimeofday(&start, NULL);
#endif

 int i1,j1,k1;
 int i,j,cnt;
 printf("@@START ITERATION :: i1, nLevels-1=%d\n", nLevels-1);
 for (i1 = 0; i1 < nLevels-1; ++i1) {
    printf("\t ITERATION :: i1 = %d #j1 From:%d To:%d\n", i1, levelPtr[i1], levelPtr[i1+1]);

#ifndef PTHREAD
   for (j1 = levelPtr[i1]; j1 < levelPtr[i1 + 1]; ++j1) {
    map = (int*)malloc(sizeof(int)* n);
    contribs = (double*)malloc(sizeof(double)* super_max * col_max);
    xi = (int*)malloc(sizeof(int)* 2*supNo);

#ifdef TIMING1
    gettimeofday(&startin, NULL);
#endif

    printf("\t\t ITERATION :: j1 = %d #k1 From:%d To:%d\n", j1, parPtr[j1], parPtr[j1+1]);
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

     top = sw_ereach_sn(supNo,cT,rT,curCol,nxtCol,col2Sup, aTree,xi,xi+supNo);
     assert(top>=0);
     PRINT("\t\t\t ITERATION :: k1 = %d #i From:%d To:%d\n", k1, top, supNo);
     for(i = top; i < supNo; ++i){
      int lSN = xi[i];

  #if DEBUG
      if(xi[top++] != lSN)
                         printf("fail");
  #endif

#else
     PRINT("\t\t\t ITERATION :: k1 = %d #i PRUNE From:%d To:%d\n", k1, prunePtr[s-1], prunePtr[s]);
     for (i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
      int lSN = pruneSet[i];
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
        assert(!isnan(cur[col * nSupR + map[cRow]]));
        assert(!isnan(contribs[i2 * nSupRs + j2]));
        assert(isFiniteNumber(cur[col * nSupR + map[cRow]]));
        assert(isFiniteNumber(contribs[i2 * nSupRs + j2]));
        cur[col * nSupR + map[cRow]] -= contribs[i2 * nSupRs + j2];
       }
      }
      PRINT("inner step : 4, %d\n", i);
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

    free(xi);
    free(contribs);
    free(map);
   } // end of j1

#else //have PTHREAD
{
    printf("prepare pthread param\n");

    int levelBgn = levelPtr[i1];
    int levelEnd = levelPtr[i1+1];
    int levelNUM = levelEnd - levelBgn;
    int cntcnt = (levelNUM-1)/NUM_CG + 1;

    printf("levelBgn = %d, levelEnd = %d, levelNUM = %d, cntcnt = %d\n", levelBgn, levelEnd, levelNUM, cntcnt);

    for (iii = 0; iii < NUM_CG; iii++){
        params[iii]->levelBgn = levelBgn + iii * cntcnt;
        params[iii]->levelEnd = (iii+1)*cntcnt > levelNUM ? levelEnd : levelBgn+(iii+1)*cntcnt;
    }

    printf("pthread param finish\n");


    for(iii = 0; iii < NUM_CG; iii++){

        int rc = pthread_create(&pthread_handler[iii], NULL, sw_cholesky_left_par_05_pthread_inner, (void*)params[iii]);
        if (rc != 0)
            printf("ERROR ##FINISH pthread create%d, rc=%d\n", iii, rc);
        //else
        //    printf("##FINISH pthread create%d, rc=%d\n", iii, rc);
    }
    for(iii = 0; iii < NUM_CG; iii++){
        int rc = pthread_join(pthread_handler[iii], NULL);
        if (rc != 0)
            printf("ERROR ##FINISH pthread join%d, rc=%d\n", iii, rc);
    }

    printf("##FINISH pthread i1 = %d\n", i1);

}//pthread
#endif

 } //end of i1
printf("/*********END OF ALL FOR**********/\n");


#ifdef TIMING
  // record the whole running time, part 0
  gettimeofday(&end, NULL);
  elapsed_seconds = TIME(start, end);
  duration2 = elapsed_seconds;
  //std::cout<<duration2<<"; ";
  timing[0]=duration2;
#endif


#ifdef TIMING
   // record the whole running time, part 1
   gettimeofday(&start, NULL);
#endif


#if 1
 printf("@@START ITERATION :: i1, last Level\n");
#ifdef TIMING
  gettimeofday(&start, NULL);
#endif
  map = (int*)malloc(sizeof(int)* n);
  contribs = (double*)malloc(sizeof(double)* super_max * col_max);
  xi = (int*)malloc(sizeof(int)* 2*supNo);
    printf("\t ITERATION :: i1 = %d #j1 From:%d To:%d\n", i1, levelPtr[nLevels-1], levelPtr[nLevels]);
  for (j1 = levelPtr[nLevels - 1]; j1 < levelPtr[nLevels]; ++j1) {
    printf("\t\t ITERATION :: j1 = %d #k1 From:%d To:%d\n", j1, parPtr[j1], parPtr[j1+1]);
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
    top = sw_ereach_sn(supNo, cT, rT, curCol, nxtCol, col2Sup, aTree, xi, xi + supNo);
     PRINT("\t\t\t ITERATION :: k1 = %d #i From:%d To:%d\n", k1, top, supNo);
    for (i = top; i < supNo; ++i) {
     int lSN = xi[i];
#else
     PRINT("\t\t\t ITERATION :: k1 = %d #i PRUNE From:%d To:%d\n", k1, prunePtr[s-1], prunePtr[s]);
     for (i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
       int lSN = pruneSet[i];
#endif

     MKL_INT nSupRs = 0;
     int cSN = blockSet[lSN];//first col of current SN
     int cNSN = blockSet[lSN + 1];//first col of Next SN
     MKL_INT Li_ptr_cNSN = Li_ptr[cNSN];
     MKL_INT Li_ptr_cSN = Li_ptr[cSN];
     MKL_INT nSNRCur = Li_ptr_cNSN - Li_ptr_cSN;
     MKL_INT supWdts = cNSN - cSN;//The width of current src SN
     int lb = 0, ub = 0;
     int sw = 1;
     for (j = Li_ptr_cSN; j < Li_ptr_cNSN; ++j) {
      //finding the overlap between curCol and curCol+supWdt in the src col
      if (Ls[j] >= curCol && sw) {
       //src*transpose(row Ls[j])
       lb = j - Li_ptr_cSN;
       sw = 0;
      }
      if (Ls[j] < curCol + supWdt && !sw) {
       ub = j - Li_ptr_cSN;
      }
     }
     nSupRs = Li_ptr_cNSN - Li_ptr_cSN - lb;
     MKL_INT ndrow1 = ub - lb + 1;
     MKL_INT ndrow3 = nSupRs - ndrow1;
     src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
     double *srcL = &lValues[lC[cSN] + ub + 1];
#ifdef BLASTIMING
     gettimeofday(&startBlas, NULL);
#endif
      PRINT("inner step : 1, %d\n", i);

#ifdef XMATH
     UPLO = 'L';
     Tran = 'N';
     dsyrk_(&UPLO,&Tran,&ndrow1,&supWdts,&one,src,&nSNRCur,&zero,
                     contribs,&nSupRs);
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
      Trana = 'N';
      Tranb = 'C';
      dgemm_(&Trana,&Tranb,&ndrow3,&ndrow1,&supWdts,&one,srcL,&nSNRCur,
                        src,&nSNRCur,&zero,contribs+ndrow1,&nSupRs );
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
     }

      PRINT("inner step : 3, %d\n", i);

     //copying contrib to L
     int i2, j2;
     for (i2 = 0; i2 < ndrow1; ++i2) {//Copy contribs to L
      int col = map[Ls[Li_ptr_cSN + i2 + lb]];//col in the SN
      for (j2 = i2; j2 < nSupRs; ++j2) {
       int cRow = Ls[Li_ptr_cSN + j2 + lb];//corresponding row in SN
       //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
        assert(!isnan(cur[col * nSupR + map[cRow]]));
        assert(!isnan(contribs[i2 * nSupRs + j2]));
        assert(isFiniteNumber(cur[col * nSupR + map[cRow]]));
        assert(isFiniteNumber(contribs[i2 * nSupRs + j2]));
       cur[col * nSupR + map[cRow]] -= contribs[i2 * nSupRs + j2];
      }
     }
      PRINT("inner step : 4, %d\n", i);
#ifndef PRUNE
     } /// end of pruneSet
#else
     } /// end of none pruneSet
#endif

#ifdef BLASTIMING
    gettimeofday(&startBlas, NULL);
#endif
#ifdef XMATH
    UPLO = 'L';
    dpotrf_(&UPLO,&supWdt,cur,&nSupR,&info);
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
    SIDE = 'R';
    UPLO = 'L';
    Trana = 'C';
    DIAG = 'N';
    dtrsm_(&SIDE, &UPLO, &Trana, &DIAG, &rowNo, &supWdt,&one,
                       cur,&nSupR,&cur[supWdt],&nSupR);
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

   }


    free(xi);
    free(contribs);
    free(map);
  } // end of j1

#endif // "#if 1"



#ifdef PTHREAD  
    for(iii = 0; iii < NUM_CG; iii++){
        free(params[iii]);
    }
#endif


#ifdef TIMING
   // record the whole running time, part 1
  gettimeofday(&end, NULL);
  elapsed_seconds = TIME(start, end);
  duration2 += elapsed_seconds;
  timing[1] = duration2;
  printf("###TOTAL TIME INSIDE numeric is %f\n\n", duration2);
#endif

#ifdef BLASTIMING
    printf("###TOTAL BLAS TIME is %f\n", blasTime);
    printf("###SYRK TIME is %f\n", syrkTime);
    printf("###GEMM TIME is %f\n", gemmTime);
    printf("###POTRF TIME is %f\n", potrfTime);
    printf("###TRSM TIME is %f\n", trsmTime);
#endif

 return 1; // return true
}

