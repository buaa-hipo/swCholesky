#ifndef SW_TRIPLE_QUEUE
#define SW_TRIPLE_QUEUE

typedef struct fusionQueue {
    int isLast;
    int ndrow1, ndrow3, supWdts;
    //int nSupRs; //nSupRs = ndrow1+ndrow3
    int nSNRCur;
    int Li_ptr_cSN;
    int lb, ub;
    double *src, *srcL;
    double *contribs;
}fusionQueue;

#ifdef __cplusplus
extern "C" {
#endif

int chooseQueue(const int ndrow1, const int ndrow3, const int supWdts);

int addTask(fusionQueue Q[64], int *cntQ, \
            const int ndrow1, const int ndrow3, const int nSupRs, \
            const int supWdts, const int nSNRCur, \
            const int Li_ptr_cSN, const int lb, const int ub, \
            double* src, double* srcL, double* contribs);

int execQueue(const fusionQueue Q[64], fusionDgemmDsyrkParam S[64], int *num, 
                const int nSupR, const int *map, const int *Ls, double *cur);

#ifdef __cplusplus
}
#endif

#endif
