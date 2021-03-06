#include <stdio.h>
#include <stdlib.h>
#include <athread.h>
#include <sys/time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "./include/common_slave.h"
#include "cblas.h"

//#define DEBUG_VERBOSE
//#define CHECK_PAD_RESULT
#define SELECT_BLK

extern void SLAVE_FUN(dgemm_dma_trans)();
extern void SLAVE_FUN(dgemm_dma_trans_alpham1_beta1)();
extern void SLAVE_FUN(copy_border_double64)();
extern void SLAVE_FUN(copy_border_back_double64)();

// M is the continous dimension
static int check_equal_val(double *a, int RM, int RK, int Me, int Ke, int Ms, int Ks, int lda)
{
  printf("%d %d -- %d %d\n", RM, RK, Me, Ke);
  int i,j;
  int cnt = 0;

  //printf("(127 128) is %lf\n", a[127*Me + 128]);
  for (j = 0; j < Ke; j++){
    for (i = 0; i < Me; i++){
      if (i >= 0 && i < Ms && j >= 0 && j < Ks) continue;
      if (i >= Ms && i < RM && j >= 0 && j < RK) {
        if (a[j*Me + i] != 1. + (j*lda+i)%100*0.01){
          printf("a(%d %d) = %lf, should=%lf\n", j, i, a[j*Me + i], 1. + (j*lda+i)%100*0.01);
          cnt++;
        }
        continue;
      }
      if (i >= 0 && i < Ms && j >= Ks && j < RK) {
        if (a[j*Me + i] != 1. + (j*lda+i)%100*0.01){
          printf("a(%d %d) = %lf, should=%lf\n", j, i, a[j*Me + i], 1. + (j*lda+i)%100*0.01);
          cnt++;
        }
        continue;
      }      

        if (a[j*Me + i] != 0.){
          printf("a(%d %d) = %lf, should=0\n", j, i, a[j*Me + i]);
          cnt++;
        }
    }
    //a[127*Me + 128] = 1;
  }
  return cnt;

}

static int check_value(const double *b, int K, int N, int ldb){
  int i, j;
  int cnt = 0;
  for (i = 0; i < K; i++){
    for (j = 0; j < N; j++){
      int off = i*ldb + j;
      if (b[off] != 1.+off%100*0.01){
          printf("b(%d %d) = %lf, should=%lf\n", i, j, b[off], 1. + (i*ldb+j)%100*0.01);
          cnt++;
      }
    }
  }
  return cnt;
}

static double estimite_compute_time(int blkM, int blkN, int blkK, int M, int N, int K) {
  //performance model achieves an estimited performance
  int bsizeN = blkN/8*sizeof(double);
  int bsizeM = blkM/8*sizeof(double);

  double a = 9.55371467e-09;
  double b = 4.80294349e-10;
  double c = 3.85210279e-11;
  double d = 1.36105221e-05;
  double T_compute = (a*blkN + b*blkM*blkN + c*blkM*blkK*blkN + d)/10
    *M/blkM*K/blkK*N/blkN;
  return T_compute;
}

/*******
 * 2019.01.06
 * A sunway-based BLAS implementation
 * suport 2 cases: Notrans A, Notrans B, Trans A, Notrans B
 * otherwise you should explicitly transpose matrix for better performance
 * Use performance model to estimite best block size
 * light-weight padding technique is applied
 *
 * Debug, I use no padding status for small case
 * *****/
void sw_cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
                 const int K, const double alpha, const double *A,
                 const int lda, const double *B, const int ldb,
                 const double beta, double *C, const int ldc)
{
  struct timeval t0, t1, t2;
  ConvData* params = (ConvData*)malloc(sizeof(ConvData));
  CopyData* cd = (CopyData*)malloc(sizeof(CopyData));


  assert(Order == CblasRowMajor);
  assert(TransB == CblasNoTrans);
  assert(TransA == CblasTrans);
  assert( (alpha == 1.0 && beta == 0.) || \
          (alpha == -1.0 && beta == 1.0) );

  assert(TransB == CblasNoTrans );
    if(TransA == CblasNoTrans)
      assert(lda >= K && ldb >= N && ldc >= N);
    else if(TransA == CblasTrans)
      assert(lda >= M && ldb >= N && ldc >= N);

    //TODO 
    //serach for best align_size and blk_size 
    //align_size is equal to blk size
    //bsizeN range from 64B ~ 2KB 
    //bsizeM range from 16B ~ 2KB
    //MBW_map should be indexed as 16Bx

#ifdef DEBUG_VERBOSE
    printf("M %d K %d N %d\n", M, K, N);
#endif
    int blkN, blkM, blkK;
    int real_blkN, real_blkM, real_blkK;
    real_blkN = 128; real_blkM = 32; real_blkK = 32;
#ifdef SELECT_BLK
    double MBW_map[]={3362.3000000000002, 6342.6000000000004, 9091.3999999999996, 11966.799999999999, 14464.4, 10109.4, 10826.799999999999, 13355.9, 14225.6, 16268.0, 17285.200000000001, 19322.400000000001, 20039.099999999999, 8748.6000000000004, 16397.0, 17568.099999999999, 18846.599999999999, 19078.799999999999, 17884.799999999999, 21040.299999999999, 21277.799999999999, 18181.299999999999, 18960.400000000001, 19724.799999999999, 20330.599999999999, 21263.700000000001, 21535.799999999999, 11486.1, 22908.099999999999, 19666.900000000001, 20302.700000000001, 21102.5, 21682.700000000001, 21875.700000000001, 22555.200000000001, 23501.799999999999, 21774.299999999999, 20105.700000000001, 21358.700000000001, 21932.099999999999, 21482.5, 19173.299999999999, 22579.200000000001, 23836.799999999999, 23775.400000000001, 21602.5, 21919.700000000001, 22429.5, 22826.400000000001, 23273.799999999999, 23630.599999999999, 24175.700000000001, 24429.099999999999, 22369.0, 21537.200000000001, 20850.400000000001, 21515.099999999999, 23762.900000000001, 23600.599999999999, 24484.400000000001, 24604.900000000001, 22800.0, 22921.599999999999, 23484.599999999999, 3390.3000000000002, 6033.3000000000002, 9180.6000000000004, 11876.4, 13766.799999999999, 10086.6, 10522.6, 13365.1, 14133.6, 16451.5, 17313.799999999999, 19445.799999999999, 19996.900000000001, 10185.9, 16442.299999999999, 17509.299999999999, 18274.200000000001, 19345.099999999999, 19300.099999999999, 20965.700000000001, 20429.400000000001, 18137.299999999999, 18710.0, 19884.799999999999, 20117.299999999999, 21051.299999999999, 21009.099999999999, 12401.799999999999, 22579.200000000001, 19770.599999999999, 20368.599999999999, 21030.200000000001, 21637.900000000001, 22160.900000000001, 22898.099999999999, 23052.599999999999, 22213.900000000001, 18106.900000000001, 21088.0, 21992.099999999999, 22231.400000000001, 18906.799999999999, 22901.299999999999, 23413.0, 23620.400000000001, 21629.400000000001, 21916.799999999999, 22201.5, 22933.200000000001, 23258.599999999999, 23628.900000000001, 24179.299999999999, 24568.700000000001, 22159.0, 22024.0, 21703.700000000001, 22000.900000000001, 23542.5, 23898.5, 24653.200000000001, 24752.200000000001, 22704.700000000001, 22908.099999999999, 23676.799999999999, 23607.799999999999, 24061.099999999999, 24265.099999999999, 24463.0, 24885.200000000001, 22663.700000000001, 23200.099999999999, 23771.5, 23822.200000000001, 23878.200000000001, 24192.599999999999, 24364.599999999999, 24847.700000000001, 23300.099999999999, 23664.5, 23730.5, 24296.799999999999, 24190.299999999999, 24188.599999999999, 24158.900000000001, 24806.599999999999, 23510.599999999999, 23840.200000000001, 24218.299999999999, 24235.200000000001, 24448.099999999999, 24828.099999999999, 25147.5, 25309.599999999999, 23574.0, 23646.700000000001, 23983.700000000001, 24372.200000000001, 24844.0, 24924.700000000001, 25109.700000000001, 25378.700000000001, 23925.400000000001, 24214.799999999999, 24551.299999999999, 24587.900000000001, 24818.299999999999, 25051.5, 25338.299999999999, 25317.400000000001, 24099.799999999999, 24158.900000000001, 24436.799999999999, 23581.200000000001, 24539.900000000001, 24685.900000000001, 24942.099999999999, 25381.200000000001, 24162.900000000001, 24270.900000000001, 24740.599999999999, 24684.099999999999, 24816.400000000001, 25079.599999999999, 25461.5, 25446.700000000001, 23707.400000000001, 23754.0, 24641.700000000001, 24930.400000000001, 25199.900000000001, 25365.099999999999, 25548.299999999999, 25656.099999999999, 24629.700000000001, 24703.0, 24873.5, 24806.599999999999, 25210.400000000001, 24867.400000000001, 24863.200000000001, 25700.5, 24606.799999999999, 24904.5, 24893.200000000001, 25083.400000000001, 25196.0, 25277.099999999999, 25601.799999999999, 24222.299999999999, 24871.599999999999, 24890.799999999999, 24978.0, 24897.0, 24061.099999999999, 25350.400000000001, 25539.799999999999, 25816.200000000001, 24777.299999999999, 25125.5, 24796.799999999999, 25261.099999999999, 25498.400000000001, 25510.700000000001, 25738.099999999999, 25643.599999999999, 24911.099999999999, 24659.200000000001, 25003.599999999999, 25063.900000000001, 25409.5, 25566.599999999999, 26173.099999999999, 25681.5, 24767.5, 25033.900000000001, 25204.200000000001, 25141.799999999999, 25419.799999999999, 25579.5, 25675.5, 25702.5, 25407.5, 24873.5, 25331.900000000001, 25430.599999999999, 24968.099999999999, 24645.400000000001, 25573.0, 25950.200000000001, 25004.099999999999, 25081.5, 25333.900000000001, 25364.599999999999, 25558.599999999999, 25914.599999999999, 25814.200000000001, 25803.599999999999, 25248.599999999999, 25129.799999999999, 25356.799999999999, 25370.900000000001, 25616.700000000001, 25525.0, 25742.700000000001, 25208.0, 25129.299999999999, 25232.599999999999, 25451.200000000001, 25446.700000000001, 25710.5, 25742.700000000001, 25948.700000000001, 25412.0, 25397.299999999999, 25393.400000000001, 25510.700000000001, 25635.599999999999, 25446.700000000001, 25925.299999999999, 25914.099999999999, 25230.700000000001, 25107.299999999999, 23748.400000000001, 25196.0, 24577.900000000001, 25139.400000000001, 25759.200000000001, 25725.599999999999, 25258.700000000001, 25432.0, 25490.0, 25510.700000000001, 25539.400000000001, 25746.700000000001, 25806.099999999999, 25706.5, 25155.700000000001, 25206.099999999999, 25381.200000000001, 25307.200000000001, 25506.299999999999, 25769.799999999999, 25683.5, 25912.0, 24691.5, 25299.400000000001, 25432.5, 25417.799999999999, 25702.5, 25772.299999999999, 24928.5, 26051.700000000001, 25049.599999999999, 25459.5, 25444.799999999999};
    double est_best_time = 1000000;
    for (blkN = 128; blkN <= N && blkN <= 2048; blkN += 128)
      for (blkM = 32; blkM <= M && blkM <= 2048; blkM += 32)
        for (blkK = 32; blkK <= K && blkK <= 2048; blkK += 32) 
    {
          int ldm_use = sizeof(double)*(blkN*blkK*2 + 
              blkK*blkM*2 + blkN*blkM)/64;
          if (ldm_use < 60*1024) {
            //printf("blkN = %d, blkM = %d, blkK = %d\n", blkN, blkM, blkK);
    
            int tmpN = (N+blkN-1)/blkN*blkN, tmpM = (M+blkM-1)/blkM*blkM, tmpK = (K+blkK-1)/blkK*blkK;

            double gflop = 2.0*N*K*M;

            int bsizeN = blkN/8*sizeof(double);
            int bsizeM = blkM/8*sizeof(double);
            double T_dma = tmpN/blkN*tmpM/blkM*tmpK/blkK*(1.0*blkN*blkK*sizeof(double)/1e6/MBW_map[bsizeN/16 - 1] + 
              1.0*blkM*blkK*sizeof(double)/1e6/MBW_map[bsizeN/16 - 1]) +
              1.0*tmpN/blkN*tmpM/blkM*blkM*blkN*sizeof(double)/1e6/MBW_map[bsizeN/16 - 1];

            double T_init_dma = (1.0*blkN*blkK*sizeof(double)/1e6/MBW_map[bsizeN/16-1] + 
              1.0*blkM*blkK*sizeof(double)/1e6/MBW_map[bsizeM/16-1]);
            double T_compute = estimite_compute_time(blkM, blkN, blkK, M, N, K);
            double est_time = MAX(T_compute, T_dma) + T_init_dma;

            if(est_time < est_best_time) {
              est_best_time = est_time;

              real_blkN = blkN;
              real_blkM = blkM;
              real_blkK = blkK;
            }
          }
    }
#endif

    blkM = real_blkM;
    blkN = real_blkN;
    blkK = real_blkK;
    // bound begin
    int Ms = (N/blkN)*blkN, Ns = (M/blkM)*blkM, Ks = (K/blkK)*blkK;
    // bound end
    int Me = (N+blkN-1)/blkN*blkN, Ne = (M+blkM-1)/blkM*blkM, Ke = (K+blkK-1)/blkK*blkK;
    int RM = N, RN = M, RK = K;

    //double* Ap = (double*)malloc(sizeof(double)*K*M*4);
    //double* Bp = (double*)malloc(sizeof(double)*K*N*4);
    //double* Cp = (double*)malloc(sizeof(double)*M*N*4);
    // CAUTION :: M N are changed
    double* Ap_ = (double*)malloc(sizeof(double)*Ke*Ne+32);
    double* Bp_ = (double*)malloc(sizeof(double)*Ke*Me+32);
    double* Cp_ = (double*)malloc(sizeof(double)*Ne*Me+32);
    double* Ap = ALIGNED(Ap_);
    double* Bp = ALIGNED(Bp_);
    double* Cp = ALIGNED(Cp_);
    params->input = B;
    params->weight = A;
    params->output = C;
    params->inputp = Bp;
    params->weightp = Ap;
    params->outputp = Cp;
    params->K = K;
    params->N = M;
    params->M = N;

    params->ldi = ldb;
    params->ldw = lda;
    params->ldo = ldc;

    params->blkN = real_blkM;
    params->blkM = real_blkN;
    params->blkK = real_blkK;



#ifdef DEBUG_VERBOSE
    printf("real_blkN %d real_blkM %d real_blkK %d\n",
        real_blkN, real_blkM, real_blkK);
    printf("RM=%d, Ms=%d, Me=%d\n", RM,Ms,Me);
    printf("RN=%d, Ns=%d, Ne=%d\n", RN,Ns,Ne);
    printf("RK=%d, Ks=%d, Ke=%d\n", RK,Ks,Ke);
    gettimeofday(&t1, NULL);
#endif

#ifdef CHECK_PAD_RESULT    
    printf("#before padding  ERROR checking B, %d\n", check_value(B, RK, RM, ldb) );
    printf("#before padding  ERROR checking A, %d\n", check_value(A, RK, RN, lda) );
#endif    

    cd->src = B;
    cd->dst = Bp;
    cd->Ms = Ms;
    cd->Ns = Ks;
    cd->M = RM;
    cd->N = RK;
    cd->Me = Me;
    cd->Ne = Ke;
    cd->blkM = blkN; //distinguish
    cd->blkN = blkK;
    cd->trans = 0;
    cd->ldx = ldb;
    athread_spawn(copy_border_double64, cd);
    athread_join();
#ifdef CHECK_PAD_RESULT    
    printf("RM=%d Ms=%d Me=%d \n", RM, Ms, Me);
    printf("RK=%d Ks=%d Ke=%d \n", RK, Ks, Ke);
    printf("ldb = %d\n", ldb);

    printf("ERROR checking B, %d\n", check_value(B, RK, RM, ldb) );
    printf("padding B\n");
    printf("ERROR padding B, %d\n", check_equal_val(Bp, RM, RK, Me, Ke, Ms, Ks, ldb) );
    printf("end padding B\n");
    //exit(0);
#endif    

#ifdef DEBUG_VERBOSE
    gettimeofday(&t2, NULL);
    double padb_time = TIME(t1,t2);
    gettimeofday(&t1, NULL);
#endif
    cd->src = A;
    cd->dst = Ap;
    cd->Ms = Ks;
    cd->Ns = Ns;
    cd->M = RK;
    cd->N = RN;
    cd->Me = Ke;
    cd->Ne = Ne;
    cd->blkM = blkK;
    cd->blkN = blkM; // distinguish
    cd->ldx = lda;
    if(TransA == CblasTrans)
      cd->trans = 1;
    else
      cd->trans = 0;
    athread_spawn(copy_border_double64, cd);
    athread_join();
#ifdef CHECK_PAD_RESULT    
    printf("ERROR checking A, %d\n", check_value(A, RK, RN, lda) );
    printf("padding A\n");
    //printf("ERROR padding A, %d\n", check_equal_val(Ap, RK, RN, Ke, Ne, Ks, Ns) );
    printf("ERROR padding A, %d\n", check_equal_val(Ap, RN, RK, Ne, Ke, Ns, Ks, lda) );
    printf("end padding A\n");
    //exit(0);
#endif

#ifdef DEBUG_VERBOSE
    gettimeofday(&t2, NULL);
    double pada_time = TIME(t1,t2);
    gettimeofday(&t1, NULL);
#endif

  if(beta==1.0){
    cd->src = C;
    cd->dst = Cp;
    cd->Ms = Ms;
    cd->Ns = Ns;
    cd->M = RM;
    cd->N = RN;
    cd->Me = Me;
    cd->Ne = Ne;
    cd->blkM = blkN; // distinguish
    cd->blkN = blkM; // distinguish
    cd->trans = 0;
    cd->ldx = ldc;
    athread_spawn(copy_border_double64, cd);
    athread_join();
#ifdef CHECK_PAD_RESULT    
    printf("ERROR checking C, %d\n", check_value(C, RN, RM, ldc) );
    printf("padding C\n");
    printf("ERROR padding C, %d\n", check_equal_val(Cp, RM, RN, Me, Ne, Ms, Ns, ldc) );
    printf("end padding C\n");
    //exit(0);
#endif
  }

#ifdef DEBUG_VERBOSE
    gettimeofday(&t2, NULL);
    double padc_time = TIME(t1,t2);
#endif

    gettimeofday(&t1, NULL);
    if(TransA == CblasNoTrans){
      //assert(lda >= K && ldb >= N && ldc >= N);
      //athread_spawn(dgemm_dma, params); }
      assert (1 && "not implemented");
    }
    else if(TransA == CblasTrans){
      //assert(lda >= M && ldb >= N && ldc >= N);
      if (alpha == 1.0 && beta == 0.)
        athread_spawn(dgemm_dma_trans, params); 
      else if (alpha == -1.0 && beta == 1.0)
        athread_spawn(dgemm_dma_trans_alpham1_beta1, params); 
    }
    athread_join();
#ifdef DEBUG_VERBOSE
    printf("FINISH calc\n");
    gettimeofday(&t2, NULL);
    double comput_time = TIME(t1,t2);
    gettimeofday(&t1, NULL);
#endif
    /*for(j = 0;j < RN; j ++)
      for(i = 0;i < RM; i ++)
        if(j >= Ns || i >= Ms)
          C[i + j*RM] = Cp[i + j*Me];*/
    cd->src = Cp;
    cd->dst = C;
    cd->Ms = Ms;
    cd->Ns = Ns;
    cd->M = RM;
    cd->N = RN;
    cd->Me = Me;
    cd->Ne = Ne;
    cd->blkM = blkN; // distinguish
    cd->blkN = blkM; // distinguish
    cd->trans = 0;
    cd->ldx = ldc;

#ifdef DEBUG_VERBOSE
    printf("[INFO] Padding start and end: Ms %d Ns %d M %d N %d Me %d Ne %d\n", cd->Ms, cd->Ns, cd->M, cd->N, cd->Me, cd->Ne);
#endif
    athread_spawn(copy_border_back_double64, cd);
    athread_join();
    /*if(Ne > Ns)
    {
      for(j = Ns; j < RN; j ++)
        memcpy(&C[j*RM], &Cp[j*Me], sizeof(float)*RM);
    }
    if(Me > Ms)
    {
      for(j = 0; j < Ns; j ++)
      {
        memcpy(&C[j*RM+Ms], &Cp[j*Me+Ms], sizeof(float)*(RM-Ms));
      }
    }*/

#ifdef DEBUG_VERBOSE
    gettimeofday(&t2, NULL);
    double depadc_time = TIME(t1,t2);
    double total_time = pada_time + padb_time + padc_time + comput_time + depadc_time;
    printf("[INFO] %lf Gflops comput_time %f pada_time %f padb_time %f padc_time %f depadc_time %f pad_percent %f\n", 
        2.0/total_time*M*N*K/1e9, comput_time, pada_time, padb_time, padc_time, depadc_time, 1-comput_time/total_time);
#endif

  free(params);
  free(cd);
  free(Ap_);
  free(Bp_);
  free(Cp_);

  return;
}


