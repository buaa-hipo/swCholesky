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
#include "sw_LSparsity.h"

#include "sw_parallel_PB_Cholesky_wavefront.h"


//#define ANALYZE
//#define ANALYZE_THEORY
#undef DEBUG
//#undef PRUNE

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
 //int numThread = atoi(argv[2]);
 //int chunk = atoi(argv[3]);
 //int costParam = atoi(argv[4]);//Inner parts
 //int levelParam = atoi(argv[5]);// level distance
 //int blasThreads = atoi(argv[6]);
 //int finalSeqNode = atoi(argv[7]);

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
 int nrelax[3] = {4,16,0};//TODO
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

 int *waveFrontPtr = new int[L->nsuper];
 int *waveFrontSet = new int[L->nsuper];
 for (int j = 0; j < L->nsuper; ++j) {
  waveFrontPtr[j]=0;
  waveFrontSet[j]=0;
 }
 nLevels = getLevelSet(L->nsuper,L->sParent,waveFrontPtr,waveFrontSet);

 printf("nLevels = %d\n", nLevels);

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

 li_ptr = L->i_ptr;
 colL = L->p;
 valL = new double[L->xsize]();
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

 gettimeofday(&start, 0L);

 cout << "START cholesky" << endl;
  printf("n = %d\n", n);
  printf("nsuper = %d\n", L->nsuper);
  printf("nLevels = %d\n", nLevels);
  printf("chunk = %d, threads = %d\n", chunk, numThread);
  printf("super_max=%d, col_max=%d\n", maxSupWid+1,maxCol+1);
  int retval = 0;
  retval = sw_cholesky_left_par_waveFront(n,A2->p,A2->i,A2->x,
                              L->p,L->s,L->i_ptr,valL,
                       L->super,L->nsuper, timingChol,
#ifndef PRUNE
                       L->sParent,A1->p, A1->i, L->col2Sup,
#else
    prunePtr,pruneSet,
#endif
                       nLevels,waveFrontPtr,waveFrontSet,
                       chunk, numThread,maxSupWid+1,maxCol+1);

 if(!retval){
    cout << "##return " << retval << endl;
    printf("###ERROR IN STEP %d\n", iterNo);
     return -1;
 }



  gettimeofday(&end, 0L);
  elapsed_seconds = TIME(start, end);
  duration2 += elapsed_seconds;
 }

 printf("@COLLECT::Annz=%d\n", NNZ);
 printf("@COLLECT::An=%d\n", n);
 printf("@COLLECT::Lnnz=%d\n", L->xsize);
 printf("@COLLECT::Time=%lf\n", duration2/iterNo);
 cout << "FINISH cholesky wavefront, time is " << duration2 << endl;



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
 delete []waveFrontPtr;
 delete []waveFrontSet;


// delete []col2sup;
#ifdef PRUNE
 delete []prunePtr; delete []pruneSet;
#endif
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
 //delete []colL;
 //delete []li_ptr;
 delete []timingChol;

 delete []inPerm;

 return 0;
}
