#include "../common/def.h"

/* find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k)) */
extern "C" 
int sw_ereach_sn (int n, int *Ap, int *Ai, int col1,int col2, int *col2sup,
               const int *parent, int *s, int *w)
{
 int i, p, len, top;
 if (!Ap || !Ai || !parent || !s || !w) return (-1) ;   /* check inputs */
 top = n;
 for (int k = col1; k < col2; ++k) {
  ASSERT(col2sup[k] < n);
  if(k==col1)
  CS_MARK (w, col2sup[k]) ;                /* mark node k as visited */
  for (p = Ap [k] ; p < Ap [k+1] ; p++){
   i = col2sup[Ai [p]] ;                /* A(i,k) is nonzero block */
   ASSERT(i < n);
   if (Ai [p] > k)
    continue ;       /* only use upper triangular part of A */
   //if(col2sup[i] == col2sup[Ai[p-1]]) continue; // from the same supernode
   for (len = 0 ; !CS_MARKED (w,i) ; i = parent [i]) /* traverse up etree*/
   {
    s [len++] = i ;         /* L(k,i) is nonzero */
    ASSERT(len < n);
    CS_MARK (w, i) ;        /* mark i as visited */
    ASSERT(i < n);
   }
   while (len > 0) s [--top] = s [--len] ; /* push path onto stack */
  }
 }
 for (p = top ; p < n ; p++) CS_MARK (w, s [p]) ;    /* unmark all nodes */
 //for (int k = col1; k < col2; ++k) {
 CS_MARK (w, col2sup[col1]);                /* unmark node k */
 //}
 return (top) ;                  /* s [top..n-1] contains pattern of L(k,:)*/
}


