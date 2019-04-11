#ifndef SWBLAS_DEF
#define SWBLAS_DEF

//#define CHOLESKY_BLAS

#define MAX_SIZE_DOUBLE 7168
#define BLAS_OFFSET(m,n,lda) ((n)*(lda)+(m))
#define ALIGNED(addr) ((((unsigned long)(addr)>>5)+1)<<5)
//static inline int BLAS_OFFSET(m,n,lda) {return n*lda+m;}
#define MIN(x,y) ((x)>(y)?(y):(x))
#define MAX(x,y) ((x)>(y)?(x):(y))

typedef struct dsyrkParam {
    int lda, ldc, k, n;
    double *a, *c;
}dsyrkParam;

typedef struct dgemmParam {
    int lda, ldb, ldc, k, m, n;
    double *a, *b, *c;
}dgemmParam;

typedef struct dtrsmParam {
    int m, n, lda, ldb;
    double *a, *b;
}dtrsmParam;

typedef struct dtrsmPadParam {
    int m, n;
    int ldb;
    int blkM;
    double *b, *bp;    
}dtrsmPadParam;

typedef struct fusionDgemmDsyrkParam {
    int ndrow1, ndrow3, wdt;
    int lda, ldc;
    double *a, *c;
}fusionDgemmDsyrkParam;

typedef struct fusionTPParam {
    fusionDgemmDsyrkParam* Q;
    int num;
}fusionTPParam;

#endif
