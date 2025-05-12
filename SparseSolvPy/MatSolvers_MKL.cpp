

#include "MatSolversEigenMKL.hpp"

#ifdef INTEL_MKL_SOLVER_USING


#include <mkl.h>
#include <omp.h>

/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
Eigenソルバ
//=======================================================
//=======================================================

*/

/*//=======================================================
// ● MKL Pardisoで解く(非対称)
//=======================================================*/
template<typename MType, typename VType>
bool MatSolversEigenMKL::solveMLKpardisoBase(const slv_int size0, const MType& matA, VType* vecB, VType *results, int mat_mode, int num_para=1){
    const slv_int n = size0;
    /* スパース形式コピー */
    slv_int *ia = new int[n+1];
    auto row_ptr = matA.matrix->matrix.outerIndexPtr();
    slv_int total;
    slv_int* ja;
    VType* a;
    if(mat_mode == 11 || mat_mode == 13){
        total = matA.matrix->matrix.nonZeros();
    }else{
        slv_int temp = matA.matrix->matrix.nonZeros();
        total = (temp-n)/2 + n;
    }
    ja = new slv_int[total];
    auto col_ptr = matA.matrix->getColPtr();
    a = new VType[total];
    auto val_ptr = matA.matrix->getValuePtr();
    /* 非対称のとき＝そのまま */
    if(mat_mode == 11 || mat_mode == 13){
        for(slv_int i = 0; i < n; ++i){
            ia[i] = row_ptr[i];
            ia[i]++;
        }
        ia[size0] = total+1;
        for (slv_int i = 0; i < total; ++i){
            ja[i] = col_ptr[i];
            ja[i]++;
            a[i] = val_ptr[i];
        }
    /* 対称のとき＝上三角だけ取り出す */
    }else{
        slv_int temp_count=0;
        slv_int push_count=0;
        for(slv_int i = 0; i < n; ++i){
            slv_int start = row_ptr[i];
            slv_int end;
            if(i == n-1){
                end = total;
            }else{
                end = row_ptr[i+1];
            }
            bool push=false;
            for(slv_int j = start ; j < end ; j++){
                slv_int retu = col_ptr[temp_count];
                VType val = val_ptr[temp_count];                
                if(retu == i){
                    push = true;
                    ia[i] = j+1;
                }
                if(push){
                    ja[push_count] = retu+1;
                    a[push_count] = val;
                    push_count++;
                }
                temp_count++;
            }

        }
        ia[size0] = total+1;
    }

    /* 行列タイプの設定設定 */
    int mtype = mat_mode;       

    /* 右辺コピー */
    VType* b = new VType[n];
    for(slv_int i = 0; i < n; ++i){
        b[i] = vecB[i];
    }

    int nrhs = 1; 
    void *pt[64];    
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;
    double ddum;    /* Double dummy */
    int idum;       /* Integer dummy. */

    omp_set_num_threads(num_para);
    for (int i = 0; i < 64; i++ ){
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[2] = num_para;
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 0;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 1;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = 1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */    
    for (int i = 0; i < 64; i++ ){
        pt[i] = 0;
    }
    /* 処理１ */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* 処理２ */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }    
    /* 処理３ */
    phase = 33;
    iparm[7] = 2;   /* Max numbers of iterative refinement steps. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, results, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* メモリ開放１ */
    phase = -1;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;

    return true;
}

/*//=======================================================
// ● MKL Pardisoで解く(対称)
//=======================================================*/
bool MatSolversEigenMKL::solveMLKpardisoSym(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para=1){
    bool bl = solveMLKpardisoBase<SparseMat, double>(size0, matA, vecB, results, 1, num_para);
    return bl;
}

/*//=======================================================
// ● MKL Pardisoで解く(複素対称)
//=======================================================*/
bool MatSolversEigenMKL::solveMLKpardisoSym(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex *results, int num_para=1){
    bool bl = solveMLKpardisoBase<SparseMatC, dcomplex>(size0, matA, vecB, results, 3, num_para);
    return bl;
}

/*//=======================================================
// ● MKL Pardisoで解く(非対称)
//=======================================================*/
bool MatSolversEigenMKL::solveMLKpardiso(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para=1){
    bool bl = solveMLKpardisoBase<SparseMat, double>(size0, matA, vecB, results, 11, num_para);
    return bl;
}

/*//=======================================================
// ● MKL Pardisoで解く(非複素対称)
//=======================================================*/
bool MatSolversEigenMKL::solveMLKpardiso(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex *results, int num_para=1){
    bool bl = solveMLKpardisoBase<SparseMatC, dcomplex>(size0, matA, vecB, results, 13, num_para);
    return bl;
}


#ifdef AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
/*//=======================================================
// ● MKL Pardisoで解く(非対称)
//=======================================================*/
bool MatSolvers::solveMLKpardiso(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para){
    const int n = size0;
    /* スパース形式コピー */
    slv_int *ia = new int[n+1];
    auto row_ptr = matA.matrix->matrix.outerIndexPtr();
    const slv_int total = matA.matrix->matrix.nonZeros();
    slv_int* ja = new slv_int[total];
    auto col_ptr = matA.matrix->getColPtr();
    double* a = new double[total];
    auto val_ptr = matA.matrix->getValuePtr();

    for(slv_int i = 0; i < n; ++i){
        ia[i] = row_ptr[i];
        ia[i]++;
    }
    ia[size0] = total+1;
    for (slv_int i = 0; i < total; ++i){
        ja[i] = col_ptr[i];
        ja[i]++;
        a[i] = val_ptr[i];
    }

    /* 実非対称設定 */
    int mtype = 11;       
    /* 右辺コピー */
    double* b = new double[n];
    for(slv_int i = 0; i < n; ++i){
        b[i] = vecB[i];
    }

    int nrhs = 1; 
    void *pt[64];    
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;
    double ddum;    /* Double dummy */
    int idum;       /* Integer dummy. */

    omp_set_num_threads(num_para);
    for (int i = 0; i < 64; i++ ){
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[2] = num_para;
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 0;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 1;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = 1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */    
    for (int i = 0; i < 64; i++ ){
        pt[i] = 0;
    }
    /* 処理１ */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* 処理２ */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }    
    /* 処理３ */
    phase = 33;
    iparm[7] = 2;   /* Max numbers of iterative refinement steps. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, results, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* メモリ開放１ */
    phase = -1;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;

    return true;
}


/*//=======================================================
// ● MKL Pardisoで解く(複素/非対称)
//=======================================================*/
bool MatSolvers::solveMLKpardiso(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex* results, int num_para){
    const slv_int n = size0;
    /* スパース形式コピー */
    slv_int *ia = new slv_int[n+1];
    auto row_ptr = matA.matrix->matrix.outerIndexPtr();
    const slv_int total = matA.matrix->matrix.nonZeros();
    slv_int* ja = new slv_int[total];
    auto col_ptr = matA.matrix->getColPtr();
    dcomplex* a = new dcomplex[total];
    auto val_ptr = matA.matrix->getValuePtr();

    for(slv_int i = 0; i < n; ++i){
        ia[i] = row_ptr[i];
        ia[i]++;
    }
    ia[size0] = total+1;
    for (slv_int i = 0; i < total; ++i){
        ja[i] = col_ptr[i];
        ja[i]++;
        a[i] = val_ptr[i];
    }

    /* 実非対称設定 */
    int mtype = 13;       
    /* 右辺コピー */
    dcomplex* b = new dcomplex[n];
    for(slv_int i = 0; i < n; ++i){
        b[i] = vecB[i];
    }

    int nrhs = 1; 
    void *pt[64];    
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;
    double ddum;    /* Double dummy */
    int idum;       /* Integer dummy. */

    omp_set_num_threads(num_para);
    for (int i = 0; i < 64; i++ ){
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[2] = num_para;
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 0;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 1;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = 1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */    
    for (int i = 0; i < 64; i++ ){
        pt[i] = 0;
    }
    /* 処理１ */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* 処理２ */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }    
    /* 処理３ */
    phase = 33;
    iparm[7] = 2;   /* Max numbers of iterative refinement steps. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, results, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* メモリ開放１ */
    phase = -1;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;

    return true;
}
#endif

/* end of namespace */
};


#else

namespace SRLfem{


template<typename MType, typename VType>
bool MatSolversEigenMKL::solveMLKpardisoBase(const slv_int size0, const MType& matA, VType* vecB, VType *results, int mat_mode, int num_para){
    return true;
}
bool MatSolversEigenMKL::solveMLKpardisoSym(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para){
    return true;
}
bool MatSolversEigenMKL::solveMLKpardisoSym(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex *results, int num_para){
    return true;
}
bool MatSolversEigenMKL::solveMLKpardiso(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para){
    return true;
}
bool MatSolversEigenMKL::solveMLKpardiso(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex *results, int num_para){
    return true;
}


};
#endif

