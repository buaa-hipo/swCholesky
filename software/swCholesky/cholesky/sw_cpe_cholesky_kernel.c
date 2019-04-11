#include <stdio.h>
#include <stdlib.h>
#include <slave.h>
#include <simd.h>
#include <dma.h>
#include "../common/sw_cpe_def.h"
#include "sw_cpe_cholesky_kernel.h"

void sw_cpe_cholesky_kernel(cpe_param *par){
    const int my_id = athread_get_id(-1);
    volatile int get_reply=0, put_reply=0;
    cpe_param param;
    get_reply = 0;
    athread_get(BCAST_MODE, par, &param, sizeof(cpe_param), &get_reply, 0, 0, 0);
    while(get_reply != 0);

    int n           = param.n;
    int *c          = param.c;
    int *r          = param.r;
    double *values  = param.values;

    size_t *lC      = param.lC;
    double *lValues = param.lValues;
    int *Ls         = param.Ls;
    size_t *Li_ptr  = param.Li_ptr;

    int *blockSet   = param.blockSet;
    int supNo       = param.supNo;

  #ifndef PRUNE
    asssert(0 && "need to implement PRUNE in CPE");
  #else
    int *prunePtr   = param.prunePtr;
    int *pruneSet   = param.pruneSet;
  #endif

    int nLevels     = param.nLevels;
    int *levelPtr   = param.levelPtr;
    int *levelSet   = param.levelSet;

    int nPar        = param.nPar;
    int *parPtr     = param.parPtr;
    int *partition  = param.partition;

    int super_max   = param.super_max;
    int col_max     = param.col_max;
    double *nodCost = param.nodCost;

    int parBgn      = param.parBgn;
    int parEnd      = param.parEnd;


    // START calculation
    int k1, j1, i, j;
    int cnt;
    for (k1 = parBgn; k1 < parEnd; k1 += SLAVES){
     int s = partition[k1] + 1;

     int curCol = s != 0 ? blockSet[s - 1] : 0;
     int nxtCol = blockSet[s];
     int supWdt = nxtCol - curCol;
     int nSupR = Li_ptr[nxtCol] - Li_ptr[curCol];//row size of supernode
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
     //TODO SUNWAY CPE
#else
     ///printf("\t\t\t ITERATION :: k1 = %d #i PRUNE From:%d To:%d\n", k1, prunePtr[s-1], prunePtr[s]);
     for (i = prunePtr[s - 1]; i < prunePtr[s]; ++i) {
      int lSN = pruneSet[i];
#endif

      int nSupRs = 0;
      int cSN = blockSet[lSN];//first col of current SN
      int cNSN = blockSet[lSN + 1];//first col of Next SN
      size_t Li_ptr_cNSN = Li_ptr[cNSN];
      size_t Li_ptr_cSN = Li_ptr[cSN];
      int nSNRCur = Li_ptr_cNSN - Li_ptr_cSN;
      int supWdts = cNSN - cSN;//The width of current src SN
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
      int ndrow1 = ub - lb + 1;
      int ndrow3 = nSupRs - ndrow1;
      src = &lValues[lC[cSN] + lb];//first element of src supernode starting from row lb
      double *srcL = &lValues[lC[cSN] + ub + 1];
#ifdef BLASTIMING
      gettimeofday(&startBlas, NULL);
#endif
      ///printf("inner step : 1, %d\n", i);

#ifdef XMATH
      //printf("inner dsyrk_\n");
      UPLO = 'L';
      Tran = 'N';
      dsyrk_(&UPLO,&Tran,&ndrow1,&supWdts,&one,src,&nSNRCur,&zero,
                             contribs,&nSupRs);
      //printf("inner dsyrk_ end\n");
#endif
      ///printf("inner step : 2, %d\n", i);

#ifdef MYBLAS
      //TODO
#endif

      if (ndrow3 > 0) {
#ifdef XMATH
      //printf("inner dgemm_\n");
       Trana = 'N'; 
       Tranb = 'C'; 
       dgemm_(&Trana,&Tranb,&ndrow3,&ndrow1,&supWdts,&one,srcL,&nSNRCur,
                                src,&nSNRCur,&zero,contribs+ndrow1,&nSupRs );
      //printf("inner dgemm_ end\n");
#endif
      ///printf("inner step : 3, %d\n", i);

#ifdef MYBLAS
       //TODO
#endif
#ifdef BLASTIMING
       gettimeofday(&endBlas, NULL);
       elapsed_seconds = TIME(startBlas, endBlas);
       blasTime += elapsed_seconds;
#endif
      }// end of ndrow3

      //copying contrib to L
      int i2, j2;
      for (i2 = 0; i2 < ndrow1; ++i2) {//Copy contribs to L
       int col = map[Ls[Li_ptr_cSN + i2 + lb]];//col in the SN
       for (j2 = i2; j2 < nSupRs; ++j2) {
        int cRow = Ls[Li_ptr_cSN + j2 + lb];//corresponding row in SN
        //lValues[lC[curCol+col]+ map[cRow]] -= contribs[i*nSupRs+j];
        cur[col * nSupR + map[cRow]] -= contribs[i2 * nSupRs + j2];
       }
      }
      ///printf("inner step : 4, %d\n", i);
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
     UPLO = 'L';
     dpotrf_(&UPLO,&supWdt,cur,&nSupR,&info);
     //printf("inner dpotrf end\n");
     if(info!=0)
         break;
#endif
#ifdef MYBLAS
     Cholesky_col(nSupR,supWdt,cur);
#endif

     int rowNo = nSupR - supWdt;
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
     for (int i = supWdt; i < nSupR; ++i) {
                     lSolve_dense_col(nSupR,supWdt,cur,&cur[i]);
                 }//TODO
#endif
#ifdef BLASTIMING
     gettimeofday(&endBlas, NULL);
     elapsed_seconds = TIME(startBlas, endBlas);
     blasTime += elapsed_seconds;
#endif



    }

}
