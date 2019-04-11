#ifndef SW_PARALLEL_PB_CHOLESKY_WAVEFRONT_H
#define SW_PARALLEL_PB_CHOLESKY_WAVEFRONT_H

#ifndef PRUNE
#define PRUNE
#endif


typedef struct WavefrontInnerParam{

    int pruneBgn;
    int pruneEnd;
    int nSupR;
    int curCol;
    int supWdt;
    double* cur;
    int* map;
    pthread_mutex_t* lock;

    // parameter for CSC matrix A
    int n;
    int* c;
    int* r;
    double* values;

    // parameter for CSC matrix L, which is supernodal
    size_t *lC; // col pointer
    int* lR; // indicates the actual row of each row in supernode
    size_t* Li_ptr; // index for Ls
    double* lValues;

    // pointer to the first col of each supernode
    int *blockSet;
    int supNo; // number of supernode

#ifndef PRUNE
    int *aTree;
    int *cT;
    int *rT;
    int *col2Sup;
#else
    int *prunePtr;
    int *pruneSet;
#endif

    int nLevels;
    int *levelPtr;
    int *levelSet;

    int super_max;
    int col_max;

}WavefrontInnerParam;

#ifdef __cplusplus
extern "C" {
#endif

int sw_cholesky_left_par_waveFront(int n, int* c, int* r, double* values,
                                 size_t *lC, int* lR, size_t* Li_ptr, double* lValues,
                          int *blockSet, int supNo, double *timing,
#ifndef PRUNE
                                 int *aTree, int *cT, int *rT, int *col2Sup,
#else
  int *prunePtr, int *pruneSet,
#endif
                          int nLevels, int *levelPtr, int *levelSet,
                          int chunk, int threads,int super_max
                          ,int col_max);


#ifdef __cplusplus
}
#endif

#endif                          
