typedef struct sw_cpe_partition_d{
    int cSN;
    int cNSN;
    size_t Li_ptr_cSN;
    size_t Li_ptr_cNSN;
    int lb;
    int ub;

    int lC_cSN; //offset for src
    int lC_cNSN; // offset for srcL
};

typedef struct sw_cpe_partation_s{
    int curCol;
    int nxtCol;
    size_t Li_ptr_curCol;
    size_t Li_ptr_nxtCol;
    int lC_curCol;

    int dNum; // number of partition d
    sw_cpe_partition_d *dList;

};

typedef struct sw_cpe_param{

    int parBgn;
    int parEnd;
    
    int sNum; // number of partition s
    sw_cpe_partition_s *sList;

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
    double *nodCost
}cpe_param;
