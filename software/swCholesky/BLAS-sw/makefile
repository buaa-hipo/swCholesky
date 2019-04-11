#CC = sw5cc -OPT:IEEE_arith=1 -D CHOLESKY_BLAS 
CC = sw5cc -D CHOLESKY_BLAS #-D PURE 
#-D DEBUG_SINGLE
LINK = mpif90 -hybrid
#LINK = sw5g++ -mhybrid

all : test_syrk_gemm

LINKFLAGS = -lm #-lm_slave #-lxMath
BLASLIB  = ../thirdparty/lib/cblas_LINUX.a
BLASLIB += ../thirdparty/lib/libswblas_opt.a
#LAPACKLIB = ../thirdparty/lib/xMath-BLAS-SPE.a
#LAPACKLIB += ../thirdparty/lib/xMath-lapack-SPE.a
LAPACKLIB = -lxMath

MYGEMMLIB = ../swGEMM/libswdgemm.a

#LINKFLAGS = -lm -lm_slave \
			-L/usr/sw-mpp/apps/share/lib -llapack -lblas -lgfortran
FLAGS = -D PURE

DSYRK_O = dsyrk_single.o dsyrk_single_slave.o dsyrk_cluster_master.o
DGEMM_O = dgemm_single.o dgemm_single_slave.o
DTRSM_O = dtrsm_single.o dtrsm_cluster_slave.o dtrsm_cluster_master.o \
			dgemm_cluster_master.o
DPOTRF_O = dpotrf_cluster_master.o dpotrf2_cluster_master.o dpotrf2_single.o \
			dsyrk_cluster_master.o \
			dtrsm_cluster_master.o dtrsm_cluster_slave.o \
			dgemm_cluster_master.o

LIB_O = 	dpotrf_cluster_master.o dpotrf2_cluster_master.o dpotrf2_single.o \
			dtrsm_cluster_master.o dtrsm_cluster_slave.o dtrsm_cluster_slave_pad.o dtrsm_single.o \
			dsyrk_cluster_master.o dsyrk_single_slave.o dsyrk_single.o \
			dgemm_cluster_master.o dgemm_single_slave.o dgemm_single.o \
			fusion_tp_master.o fusion_dgemm_dsyrk_single_slave.o

lib: $(LIB_O)
	sw5ar cr ./libswCHOLBLAS.a $^

test_dpotrf:test_dpotrf.o
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) -lxMath $(LINK_SPC)
	sh run.sh ./$@
test_dtrsm: test_dtrsm.o #dtrsm_cluster_master.o dtrsm_cluster_slave.o dtrsm_cluster_slave_pad.o dtrsm_single.o dgemm_cluster_master.o dgemm_single_slave.o dgemm_single.o
	#$(LINK) $^ -o $@ ./libswCHOLBLAS.a $(MYGEMMLIB) -allshare -lpthread
	#$(LINK) $^ -o $@ $(LINKFLAGS) $(MYGEMMLIB) -allshare -lpthread
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) -lxMath $(LINK_SPC)
	sh run.sh ./$@
test_dgemm: test_dgemm.o
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) $(LINK_SPC)
	sh run.sh ./$@
test_dsyrk: test_dsyrk.o
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) $(LINK_SPC)
	sh run.sh ./$@
test_syrk_gemm: test_syrk_gemm.o
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) $(LINK_SPC)
	sh run.sh ./$@
test_dgemm_plot: test_dgemm_plot.o
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) $(LINK_SPC)
	sh run.sh ./$@
test_dsyrk_plot: test_dsyrk_plot.o
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) $(LINK_SPC)
	sh run.sh ./$@
test_fusion: test_fusion.o
	$(LINK) $^ -o $@ $(LINKFLAGS) ./libswCHOLBLAS.a $(BLASLIB) $(MYGEMMLIB) $(LINK_SPC)
	sh run.sh ./$@

test_fusion.o: test_fusion.c
	$(CC) -host -c $< -o $@ $(FLAGS)
fusion_dgemm_dsyrk_single_slave.o: fusion_dgemm_dsyrk_single_slave.c
	$(CC) -slave -msimd -c $< -o $@ -O
fusion_tp_master.o: fusion_tp_master.c
	$(CC) -host -c $< -o $@ -O3

test_dtrsm.o: test_dtrsm.c
	$(CC) -host -c $< -o $@ $(FLAGS)
dtrsm_single.o: dtrsm_single.c
	$(CC) -host -c $< -o $@ -O3  
dtrsm_cluster_slave.o: dtrsm_cluster_slave.c
	$(CC) -slave -msimd -c $< -o $@ -O3
dtrsm_cluster_slave_pad.o: dtrsm_cluster_slave_pad.c
	$(CC) -slave -msimd -c $< -o $@ -O3
dtrsm_cluster_master.o: dtrsm_cluster_master.c
	$(CC) -host -c $< -o $@ -O3

test_dpotrf.o: test_dpotrf.c
	$(CC) -host -c $< -o $@ $(FLAGS)
dpotrf_cluster_master.o: dpotrf_cluster_master.c
	$(CC) -host -c $< -o $@ -O3
dpotrf2_cluster_master.o: dpotrf2_cluster_master.c
	$(CC) -host -c $< -o $@ -O3
dpotrf2_single.o: dpotrf2_single.c
	$(CC) -host -c $< -o $@ -O3  

test_dgemm_plot.o: test_dgemm_plot.c
	$(CC) -host -c $< -o $@ $(FLAGS)
test_dgemm.o: test_dgemm.c
	$(CC) -host -c $< -o $@ $(FLAGS)
dgemm_single.o: dgemm_single.c
	$(CC) -host -c $< -o $@ -O3
dgemm_single_slave.o: dgemm_single_slave.c
	$(CC) -slave -msimd -c $< -o $@ -O3
dgemm_cluster_master.o: dgemm_cluster_master.c
	$(CC) -host -c $< -o $@ -O3

test_dsyrk_plot.o: test_dsyrk_plot.c
	$(CC) -host -c $< -o $@ $(FLAGS)
test_dsyrk.o: test_dsyrk.c
	$(CC) -host -c $< -o $@ $(FLAGS)
dsyrk_single.o: dsyrk_single.c
	$(CC) -host -c $< -o $@ -O3
dsyrk_single_slave.o: dsyrk_single_slave.c
	$(CC) -slave -msimd -c $< -o $@ -O3 
dsyrk_cluster_master.o: dsyrk_cluster_master.c
	$(CC) -host -c $< -o $@ -O3

test_syrk_gemm.o: test_syrk_gemm.c
	$(CC) -host -c $< -o $@ $(FLAGS)

clean:
	rm -f $(DSYRK_O) $(DGEMM_O) $(DTRSM_O) $(DPOTRF_O) test_*.o 
