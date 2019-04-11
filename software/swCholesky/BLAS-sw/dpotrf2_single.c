#include <math.h>

//void Cholesky_col(int n, int dim, double* a){
void dpotrf2_single(const int dim, double* a, const int n, int* info){
 double tmp = 0;
 int i,j,k;
 for (j = 0; j < dim; ++j) {
  for (k = 0; k < j; ++k) {
   tmp = a[k*n+j];
   for (i = j; i < dim; ++i) {
    a[j*n+i] = a[j*n+i] - a[k*n+i]*tmp;
   }
  }
  if(a[j*n+j] <= 0 || isnan(a[j*n+j])){
      *info = j;
      return;
  }
  tmp = sqrt(a[j*n+j]);
  for (k = j+1; k < dim; ++k) {
   a[j*n+k] = a[j*n+k] / tmp;
  }
  a[j*n+j] = tmp;
 }
}
