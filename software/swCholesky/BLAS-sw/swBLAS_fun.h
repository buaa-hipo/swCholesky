#ifndef SWBLAS_FUN
#define SWBLAS_FUN

void dgemm_cluster_master_pure(const int m, const int n, const int k, const double alpha, const double* a, const int lda, const double* b, const int ldb, const double beta, double* c, const int ldc);
void dgemm_cluster_master(const int m, const int n, const int k, const double alpha, const double* a, const int lda, const double* b, const int ldb, const double beta, double* c, const int ldc);
void dgemm_single(const int m, const int n, const int k, const double alpha, const double* a, const int lda, const double* b, const int ldb, const double beta, double* c, const int ldc);
void dsyrk_cluster_master_pure(const int n, const int k, const double alpha, const double* a, const int lda, const double beta, double* c, const int ldc);
void dsyrk_cluster_master(const int n, const int k, const double alpha, const double* a, const int lda, const double beta, double* c, const int ldc);
void dsyrk_single(const int n, const int k, const double alpha, const double* a, const int lda, const double beta, double* c, const int ldc);
void dtrsm_cluster_master(const int m, const int n, const double alpha, const double* a, const int lda, double* b, const int ldb);
void dtrsm_single(const int m, const int n, const double alpha, const double* a, const int lda, double* b, const int ldb);
void dpotrf_cluster_master(const int n, double* a, const int lda, int* info);
void dpotrf2_cluster_master(const int n, double* a, const int lda, int* info);
void dpotrf2_single(const int n, double* a, const int lda, int* info);

#endif
