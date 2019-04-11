#ifndef SW_PARALLEL_PB_CHOLESKY_05_H
#define SW_PARALLEL_PB_CHOLESKY_05_H
typedef struct CholeskyInnerParam{

    int levelBgn;
    int levelEnd;

    // parameter for CSC matrix A
    int n;
    int* c;
    int* r;
    double* values;

    // parameter for CSC matrix L, which is supernodal
    size_t *lC; // col pointer
    double* lValues;
    int* Ls; // indicates the actual row of each row in supernode
    size_t* Li_ptr; // index for Ls

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

    int nPar;
    int *parPtr;
    int *partition;

    int super_max;
    int col_max;
    double *nodCost;

}InnerParam;

#ifdef __cplusplus
extern "C" {
#endif

int sw_cholesky_left_par_05(int n, int* c, int* r, double* values,
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
                          ) ;


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
                          ) ;
#ifdef __cplusplus
}
#endif

#endif
