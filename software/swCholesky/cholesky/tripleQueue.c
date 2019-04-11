
#include <stdio.h>
#include <stdlib.h>
#include <athread.h>
#include <assert.h>
#include <float.h>
#include <sys/time.h>
#include "../BLAS-sw/swBLAS_def.h"
#include "../BLAS-sw/swBLAS_fun.h"
#include "tripleQueue.h"

extern int isFiniteNumber(double a);

// return 0: JIT
// return 1: TP
// return 2: DP
inline int chooseQueue(const int ndrow1, const int ndrow3, const int supWdts){
   // gemm: m64 n32 k32
   // syrk: n16 k16
    if (ndrow1 <= 32 && ndrow3 <= 32 && supWdts <= 32)
        return 0;
    else if ((ndrow1+ndrow3)*(supWdts+1) <= 7168)
        return 1;
    else 
        return 2;

}

int addTask(fusionQueue Q[64], int *cntQ, \
            const int ndrow1, const int ndrow3, const int nSupRs, \
            const int supWdts, const int nSNRCur, \
            const int Li_ptr_cSN, const int lb, const int ub, \
            double* src, double* srcL, double* contribs)
{
    assert(srcL == src + ndrow1);
    assert(nSupRs == ndrow1+ndrow3);

    #ifdef XMATH
        char UPLO = 'L';
        char Tran = 'N';
        char Trana = 'N'; 
        char Tranb = 'C'; 
    #endif
    
    double one = 1.0, zero = 0.;
    //printf("queue case %d, (%d,%d,%d)\n", chooseQueue(ndrow1,ndrow3,supWdts), ndrow1,ndrow3,supWdts);
    //printf("src= %p, contribs= %p\n", src, contribs);
    switch(chooseQueue(ndrow1,ndrow3,supWdts)) {
        case 0: {
            dsyrk_single(ndrow1,supWdts,one,src,nSNRCur,zero,contribs,nSupRs);
            if (ndrow3 > 0)
                dgemm_single(ndrow3,ndrow1,supWdts,one, \
                    srcL,nSNRCur,src,nSNRCur,zero,contribs+ndrow1,nSupRs);
            return 0;
        }
        case 2: {
        #ifdef MYBLAS
            dsyrk_cluster_master_pure(ndrow1,supWdts,one,src,nSNRCur,zero,contribs,nSupRs);
            if (ndrow3 > 0)
                dgemm_cluster_master_pure(ndrow3,ndrow1,supWdts,one, \
                    srcL,nSNRCur,src,nSNRCur,zero,contribs+ndrow1,nSupRs);
        #endif
        #ifdef XMATH
            dsyrk_(&UPLO,&Tran,&ndrow1,&supWdts,&one,src,&nSNRCur,&zero,
                             contribs,&nSupRs);
            if (ndrow3 > 0)
                dgemm_(&Trana,&Tranb,&ndrow3,&ndrow1,&supWdts,&one,srcL,&nSNRCur,
                                src,&nSNRCur,&zero,contribs+ndrow1,&nSupRs );
        #endif
            return 2;
        }
        case 1: {
            const int n = *cntQ;
            assert(n < 64 && "DP Queue full");
            Q[n].ndrow1 = ndrow1;
            Q[n].ndrow3 = ndrow3;
            Q[n].supWdts = supWdts;
            Q[n].nSNRCur = nSNRCur;
            Q[n].Li_ptr_cSN = Li_ptr_cSN;
            Q[n].lb = lb;
            Q[n].ub = ub;
            Q[n].src = src;
            Q[n].srcL = srcL;
            Q[n].contribs = (double*)malloc(sizeof(double)*nSupRs*ndrow1);
            assert(Q[n].contribs != NULL);

            (*cntQ) += 1;            

            if ((*cntQ) < 64) return 1;
            else            return -1;
        }
        default: {
            assert(0 && "chooseQueue Error");
        }
    }
}

int execQueue(const fusionQueue Q[64], fusionDgemmDsyrkParam S[64], int *cntQ, 
                const int nSupR, const int *map, const int *Ls, double *cur)
{
    int i;
    const int num = (*cntQ);

    if (num == 0) return;

    for (i = 0; i < num; i++){
        S[i].ndrow1 = Q[i].ndrow1;
        S[i].ndrow3 = Q[i].ndrow3;
        S[i].wdt = Q[i].supWdts;
        S[i].lda = Q[i].nSNRCur;
        S[i].ldc = Q[i].ndrow1 + Q[i].ndrow3;
        S[i].a = Q[i].src;
        S[i].c = Q[i].contribs;
    }

    double one = 1.0, zero = 0.;
    if (num >= 8) {
        //printf("Q in MPE is %p\n", S);
        //printf("size in MPE = %d\n", sizeof(fusionDgemmDsyrkParam) );
        fusion_tp_master(S, num);
    }
    else {
        for (i = 0; i < num; i++){
            int ndrow1 = Q[i].ndrow1;
            int ndrow3 = Q[i].ndrow3;
            int supWdts = Q[i].supWdts;
            int nSNRCur = Q[i].nSNRCur;
            int nSupRs = ndrow1 + ndrow3;
            double *src = Q[i].src;
            double *srcL = Q[i].srcL;
            double *contribs = Q[i].contribs;
            dsyrk_single(ndrow1,supWdts,one,src,nSNRCur,zero,contribs,nSupRs);
            if (ndrow3 > 0)
                dgemm_single(ndrow3,ndrow1,supWdts,one, \
                    srcL,nSNRCur,src,nSNRCur,zero,contribs+ndrow1,nSupRs);
        }
    }
    

    int i2, j2;
    for (i = 0; i < num; i++){
        int Li_ptr_cSN = Q[i].Li_ptr_cSN;
        int lb = Q[i].lb;
        int ub = Q[i].ub;
        int nSupRs = Q[i].ndrow1 + Q[i].ndrow3;
        int ndrow1 = Q[i].ndrow1;
        double* contribs = Q[i].contribs;

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

      free(contribs);
    }

    *cntQ = 0;

    return 0;
}
