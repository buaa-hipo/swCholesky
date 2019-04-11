CC = sw5cc
CPP = mpiCC -host
#LINK = mpif90 -hybrid -lstdc++
#LINK = mpiCC -hybrid
LINK = sw5g++ -mhybrid
#LINK = mpicxx -hybrid

INCLUDEFLAGS = -I./common/ \
				-I./cholesky/

OPTFLAGS = -O2 -OPT:IEEE_arith=1
FLAGS = -D GIVEN -D VERIFY -D COLDCACHE -D CHOL  -D PRUNE \
		#-D DEBUG_CHECK 
#FLAGS += -D XMATH
FLAGS += -D XMATH -D PTHREAD
#FLAGS += -D MYBLAS 
#FLAGS += -D MYBLAS -D TRIPLEQUEUE 
#FLAGS += -D XMATH -D TRIPLEQUEUE 
#FLAGS += -D MYBLAS -D PTHREAD
#FLAGS += -D MYBLAS -D TRIPLEQUEUE -D PTHREAD
MPEFLAGS = -host  -msimd $(FLAGS)
CPEFLAGS = -slave -msimd $(FLAGS)

#LINKFLAGS = -lm -allshare -lpthread 
LINKFLAGS = -lm -lpthread -L/usr/sw-mpp/apps/share/lib -llapack -lblas -lgfortran\
			-mallshare
#BLASLIB = ./thirdparty/lib/cblas_LINUX.a
#BLASLIB += ./thirdparty/lib/libswblas_opt.a
BLASLIB = -lxMath_multicore
MYLIB = ./BLAS-sw/libswCHOLBLAS.a
#MYLIB = ./BLAS-sw/libswCHOLBLAS_xmath_combined.a
#MYLIB = ./BLAS-sw/libswCHOLBLAS_DP.a
MYLIB += ./swGEMM/libswdgemm.a
CPPFLAGS = $(FLAGS)

#CHOLESKY_O = sw_parallel_PB_Cholesky_05.o sw_parallel_PB_Cholesky_05_pthread_inner.o 
CHOLESKY_O = sw_parallel_PB_Cholesky_05.o sw_parallel_PB_Cholesky_05_pthread_inner.o \
				sw_choleskyTest01.o  sw_Reach.o def.o \
				tripleQueue.o \
				#sw_parallel_PB_Cholesky_wavefront.o sw_parallel_PB_Cholesky_wavefront_pthread_inner.o


sw_cholesky: $(CHOLESKY_O)
	$(LINK) ${OPTFLAGS} $^ -o $@ $(LINKFLAGS) $(MYLIB) 
	#$(BLASLIB)  
	#$(LINK_SPC)
sw_choleskyTest01.o : ./examples/sw_choleskyTest01.cpp
	$(CPP) ${OPTFLAGS} $(CPPFLAGS) -c $< -o $@ $(INCLUDEFLAGS)
sw_choleskyTest03.o : ./examples/sw_choleskyTest03.cpp
	$(CPP) ${OPTFLAGS} $(CPPFLAGS) -c $< -o $@ $(INCLUDEFLAGS)
def.o: ./common/def.cpp
	$(CPP) ${OPTFLAGS} $(CPPFLAGS) -c $< -o $@ $(INCLUDEFLAGS)
sw_Reach.o: ./common/sw_Reach.cpp
	$(CPP) ${OPTFLAGS} $(CPPFLAGS) -c $< -o $@ $(INCLUDEFLAGS)
sw_parallel_PB_Cholesky_05_pthread_inner.o : ./cholesky/sw_parallel_PB_Cholesky_05_pthread_inner.c
	$(CC) ${OPTFLAGS} $(MPEFLAGS) -c $< -o $@
sw_parallel_PB_Cholesky_05.o : ./cholesky/sw_parallel_PB_Cholesky_05.c
	$(CC) ${OPTFLAGS} $(MPEFLAGS) -c $< -o $@
sw_parallel_PB_Cholesky_05_check.o : ./cholesky/sw_parallel_PB_Cholesky_05_check.c
	$(CC) ${OPTFLAGS} $(MPEFLAGS) -c $< -o $@
sw_parallel_PB_Cholesky_wavefront.o : ./cholesky/sw_parallel_PB_Cholesky_wavefront.c
	$(CC) ${OPTFLAGS} $(MPEFLAGS) -c $< -o $@
sw_parallel_PB_Cholesky_wavefront_pthread_inner.o : ./cholesky/sw_parallel_PB_Cholesky_wavefront_pthread_inner.c
	$(CC) ${OPTFLAGS} $(MPEFLAGS) -c $< -o $@
tripleQueue.o: ./cholesky/tripleQueue.c
	$(CC) ${OPTFLAGS} $(MPEFLAGS) -c $< -o $@

.PHONY: clean

clean:
	rm -f $(CHOLESKY_O) sw_cholesky
