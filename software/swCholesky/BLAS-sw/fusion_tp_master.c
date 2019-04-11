#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <athread.h>
#include "swBLAS_def.h"
#include "swBLAS_fun.h"

extern SLAVE_FUN(fusion_dgemm_dsyrk_single_slave)();

void fusion_tp_master(fusionDgemmDsyrkParam Q[64], const int num){
    if (num <= 0) return;
    fusionTPParam TP;
    TP.Q = Q;
    TP.num = num;
    athread_spawn(fusion_dgemm_dsyrk_single_slave, &TP);
    athread_join();

    return;
}
