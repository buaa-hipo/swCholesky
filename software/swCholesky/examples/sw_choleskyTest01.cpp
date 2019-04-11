#include <iostream>
#include <fstream>
//#include <chrono>
#include <sys/time.h>
#include <algorithm>
//#include <cholmod.h>
//#include <cholmod_function.h>
//#include <omp.h>


#include "Ordering.h"
#include "Inspection_Prune.h"
#include "Inspection_Block.h"
#include "Util.h"
#include "PB_Cholesky.h"
#include "sw_LSparsity.h"

#include "sw_parallel_PB_Cholesky_05.h"


//#define ANALYZE
//#define ANALYZE_THEORY
#undef DEBUG

//#define CHECK_MYBLAS_RESULT

#define TIME(a,b) (1.0*((b).tv_sec-(a).tv_sec)+0.000001*((b).tv_usec-(a).tv_usec))

using namespace std;


int main(int argc, char *argv[]) {


 //std::string f1 = "/home/kazem/UFDB/SymFull/cbuckle.mtx";
 //string fName1= "/home/kazem/UFDB/SymSparsity/cbuckle_sparsity_amd.dat";

 if(argc<2)
  printf("input args are missing");
 std::string f1 = argv[1];
//    string fName1 = argv[2];
 int *col, *row;
 int  *rowL;
 size_t *colL;
 double *valL;
 double *valL_check;
 double  *y, *val, *x;
 int  maxSupWid, maxCol;
 size_t n, NNZ;
 struct timeval  start, end;
 double elapsed_seconds;
 double durationSym = 0 ,duration3 = 0, duration2=0, duration1=0;
 long totalIntraContribNNZs=0,totalInterCotribNNZs=0, numOfEdgeCuts=0;
 int numberOfIntraCore=0,numberOfInterCore=0;

 cout << "START read matrix" << endl;
 if (!readMatrix(f1,n,NNZ,col,row,val)){
  cout << "ERROR read matrix" << endl;
  return -1;
 }

/* if(wasteFile.fail())
  return -1;*/
 int chunk = atoi(argv[2]);
 int costParam = atoi(argv[3]);//Inner parts
 int levelParam = atoi(argv[4]);// level distance
 int finalSeqNode = atoi(argv[5]);
 size_t *inPerm;
 cout << "##FINISH getting param" <<endl;
 if(argc>6){
  std::string orderFileName = argv[6];
  inPerm = new size_t[n]();
  ///readOrdering(orderFileName,n,inPerm);
  FILE * pFile;
   if((pFile = fopen (orderFileName.c_str(), "rb+"))==NULL){
     printf("cant open the file\n");
     exit(0);
   }
   fread(inPerm, sizeof(size_t), n, pFile);
   fclose(pFile);
   cout << "##FINISH reading perm" <<endl;
 }


 const int numThread = 1;
 const int numThread_sw = 64;
 ///omp_set_num_threads(numThread);

 // MKL_Set_Num_Threads(1);
 ///MKL_Domain_Set_Num_Threads(blasThreads,MKL_DOMAIN_BLAS);
 //cout<<"---" <<MKL_Domain_Get_Max_Threads(MKL_DOMAIN_BLAS)<<"\n";
/*    chrono::time_point<std::chrono::system_clock> start, end;
    double durationAlltogether = 0, durationPruned=0, durationSym=0,
            ordering=0, durationBlock=0, durationAllSmall=0;
    chrono::duration<double> elapsed_seconds;*/

 int factorSize=0;
 double timing[4];//for time measurement

 ifstream spFile1;
//    spFile1.open(fName1);
 int *prunePtr, *pruneSet;
 int *levelPtr = NULL, *levelSet = NULL, *parPtr = NULL,
   *partition =NULL;
 int nLevels=0, nPar=0;

 //double timingChol[4]={.0,.0,.0,.0};//for time measurement
 double *timingChol = new double[4+numThread]();//for time measurement
 double orderingTime=0;
/* int nrelax[3] = {4,16,48};//TODO
 double zrelax[3] = {0.8,0.1,0.05};*/
 int nrelax[3] = {4,16,48};//TODO
 double zrelax[3] = {0.8,0.1,0.05};
 int status=0;
 double *contribs;
 int super_max = 164; //tunig parameter for the max size of supernodes TODO: find it in analysis
 int col_max = n;
// int *col2sup=new int[n]();
 //int *blockSet;
 //contribs = new double[super_max * col_max]();
 size_t *li_ptr; //= new size_t[n+1];
 //int *map = new int[n]();
 //colL = new size_t[n + 1]();
 CSC *Amat = new CSC;
 Amat->nzmax = NNZ; Amat->ncol=Amat->nrow=n;
 Amat->stype=-1;Amat->xtype=CHOLMOD_REAL;Amat->packed=TRUE;
 Amat->p = col; Amat->i = row; Amat->x=val; Amat->nz = NULL;
 Amat->sorted = TRUE;
 //start = std::chrono::system_clock::now();
 cout << "START analyze_p2" << endl;
 gettimeofday(&start, 0L);
 BCSC *L = analyze_p2(1,Amat,NULL,NULL,nrelax,zrelax,
                      n,prunePtr,pruneSet,
                      nLevels, levelPtr,levelSet,
                      nPar, parPtr, partition,
                      costParam,levelParam,finalSeqNode,
                      status, maxSupWid, maxCol, orderingTime,
                      inPerm);
 //end = std::chrono::system_clock::now();
 //elapsed_seconds = end-start;
 //durationSym=elapsed_seconds.count();
 gettimeofday(&end, 0L);
 elapsed_seconds = TIME(start, end);
 durationSym = elapsed_seconds;
 cout << "FINISH analyze_p2" << endl;
#if 0
 cout<<"\n";
    for (int j = 0; j < L->nsuper; ++j) {
        int colLength = L->pi[j+1]-L->pi[j];
        int supWid = L->super[j+1]-L->super[j];
        cout<<"Supernode: "<<j<<" Len and Wid are: "
            <<colLength<<","<<supWid<<"\n";
        for (int i = L->pi[j]; i < L->pi[j+1]; ++i) {
            cout<<L->s[i]<<",";
        }
        cout<<"\n";
    }
#endif
 //Some conversion for sympiler

 int colLength=0;

/* for (int j = 0; j < L->nsuper; ++j) {
  int curCol = L->super[j];
  int nxtCol = L->super[j+1];
  colLength = L->pi[j+1]-L->pi[j];
  //int supWid = L->super[j+1]-L->super[j];
  for (int i = curCol+1; i < nxtCol+1; ++i) {
   li_ptr[i-1] = L->pi[j];
   colL[i]= colL[curCol] + (i-curCol)*colLength;
  }
 }
 li_ptr[n] = L->pi[L->nsuper];
 colL[n] = colL[n-1] + colLength;*/
 li_ptr = L->i_ptr;
 colL = L->p;
 valL = new double[L->xsize]();
 valL_check = new double[L->xsize]();
 //delete []L->x;
 //delete []L->px;
 //delete []L->p;
 delete []L->pi;
 delete []L->i;
 delete []L->ColCount;
#if 0

 for (int j = 0; j < 10; ++j) {
        int curCol = L->super[j];
        int nxtCol = L->super[j+1];
        for (int k = li_ptr[curCol]; k < li_ptr[nxtCol]; ++k) {
            cout<<L->s[k]<<",";

        }
        cout<<"\n";
    }
#endif

 CSC *A1 = ptranspose(Amat,2,L->Perm,NULL,0,status);
 CSC *A2 = ptranspose(A1,2,NULL,NULL,0,status);
#ifndef VERIFY
 delete []L->Perm;
#endif
// int *map = new int[n]();
 //enableColdCache(1200,wasteFile);
 int iterNo=1;
 for (int k = 0; k < iterNo; ++k) {

  for (int i = 0; i < L->xsize; ++i) {
   valL[i]=0.0;
  }
  for (int i = 0; i < numThread + 4; ++i) {
   timingChol[i]=0.0;
  }

const int NN = 1;
double time = 0;
for(int ii = 0; ii < NN; ii++){

  bool retval=false;
 cout << "START cholesky" << endl;
 gettimeofday(&start, 0L);

  retval=sw_cholesky_left_par_05(n,A2->p,A2->i,A2->x,L->p,L->s,L->i_ptr,valL,
                       L->super,L->nsuper, timingChol,
#ifndef PRUNE
                       L->sParent,A1->p, A1->i, L->col2Sup,
#else
                       prunePtr,pruneSet,
#endif
                       nLevels,levelPtr,levelSet,
                       nPar, parPtr, partition,
                       chunk, numThread, maxSupWid+1,maxCol+1, NULL);


 gettimeofday(&end, 0L);
 elapsed_seconds = TIME(start, end);
 duration2 += elapsed_seconds;
 time += elapsed_seconds;
 printf("@NO::Time=%f\n", TIME(start,end));
 if(!retval){
    cout << "##return " << retval << endl;
    printf("###ERROR IN STEP %d\n", ii);
     return -1;
 }

 for(int i = 0; i < L->xsize; i++)
    valL[i] = 0.0;
}
 printf("@COLLECT::Annz=%d\n", NNZ);
 printf("@COLLECT::An=%d\n", n);
 printf("@COLLECT::Lnnz=%d\n", L->xsize);
 printf("@COLLECT::Time=%lf\n", time/NN);
 cout << "FINISH cholesky, time is " << duration2 << endl;

#ifdef CHECK_MYBLAS_RESULT
  bool retval_check=false;
 cout << "### CHECKING :: START cholesky check" << endl;
 gettimeofday(&start, 0L);

  retval_check=sw_cholesky_left_par_05_check(n,A2->p,A2->i,A2->x,L->p,L->s,L->i_ptr,valL_check,
                       L->super,L->nsuper, timingChol,
#ifndef PRUNE
                       L->sParent,A1->p, A1->i, L->col2Sup,
#else
                       prunePtr,pruneSet,
#endif
                       nLevels,levelPtr,levelSet,
                       nPar, parPtr, partition,
                       chunk, numThread, maxSupWid+1,maxCol+1, NULL);


 gettimeofday(&end, 0L);
 elapsed_seconds = TIME(start, end);
 duration2 = elapsed_seconds;
 cout << "### CHECKING :: FINISH cholesky check, time is " << duration2 << endl;
 if(!retval){
    cout << "##return check" << retval_check << endl;
     return -1;

 }
#endif
 }// end of iterNo


/* std::cout<<"BLAS kernel time: \n";
 double sumBlas=0;
 for (int i = 4; i < 4 + numThread; ++i) {
  sumBlas += timeArray[mid].tArray[i];
 }
  std::cout<<sumBlas<<"\n";*/

 allocateAC(Amat,0,0,0,FALSE);
 allocateAC(A1,0,0,0,FALSE);
 allocateAC(A2,0,0,0,FALSE);
 //allocateLC(L,FALSE);
 //TODO HACK
 
 delete []L->super;
 delete []L->sParent;
 delete []L->s;
 delete []L->col2Sup;
 delete []L->p;
 delete []L->i_ptr;
 
//#if DEBUG > 0
#ifdef CHECK_MYBLAS_RESULT
cout << "\n\n######## COMPARING THE RESULT #########" << endl;
 for (int i = n-10; i < n; ++i) {
            std::cout<<i<<":\n";
            for (int m = colL[i],cnt=0; m < colL[i+1]; ++m, ++cnt) {
                if(!std::isfinite(valL[m])) {
                    if(!std::isfinite(valL_check[m])) {
                        std::cout << "Error in col (both version)"<< i;
                        return -1;
                    } else {
                        std::cout << "Error in col"<< i;
                        return -1;
                    }
                }
                if(rowL[li_ptr[i]+cnt] >= i ){
                    if (fabs(valL[m]-valL_check[m]) > 1e-3)
                        std::cout<<valL[m]<<", " << valL_check[m] <<std::endl;

                }
            }
            std::cout<<"\n";
        }
        std::cout<<"\n";
cout << "######## FINISH COMPARING #########\n\n" << endl;
#endif

// delete []col2sup;
#ifdef PRUNE
 delete []prunePtr; delete []pruneSet;
#endif

 cout << "alloc3" << endl;
 if(levelPtr!=NULL)
  delete []levelPtr;
 if(levelPtr!=NULL)
  delete []levelSet;
 if(parPtr!=NULL )
  delete []parPtr;
 if(partition!=NULL)
  delete []partition;
 //delete []contribs;
 //delete []map;
 delete []valL;
 delete []valL_check;
 //delete []colL;
 //delete []li_ptr;
 delete []timingChol;

 delete []inPerm;

 return 0;
}


