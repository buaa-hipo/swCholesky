
# swCholesky

## build

### 1.build swGEMM
```bash
    cd swCholesky/swGEMM
    make lib
```

### 2. build swCHOLBLAS
```bash
    cd swCholesky/BLAS_sw
    make lib
```

### 3. build swCholesky

#### build pthread + swCHOLBLAS version
pthread number:: modify #define NUM_CG ?? , ?? should be 1,2,3,4
```bash
    cd swCholesky
    make -f makefile.all sw_cholesky
```

#### build mpe version (not support pthread)
```bash
    module load sw/compiler/gcc530
    cd swCholesky
    make -f makefile.mpe sw_cholesky
```

#### build xMath version (not support pthread)
```bash
    cd swCholesky
    make -f makefile.xmath sw_cholesky
```

## run
1. Download sparse positive definite matrices, which are of mtx format, to folder 'mtx_path'

2. Some matrices may be too large to run on one node, because of limited main memory

3. Generate metis permutation files using a x86 machine, to folder 'perm_path'

4. Modify parameters in test.py. mtx_path = 'mtx_path' && perm_x86_path = 'perm_path' 

```bash
    cd swCholesky
    python test.py
```

## 
swCholesky is based on parsy_bench(https://github.com/cheshmi/parsy_bench).
