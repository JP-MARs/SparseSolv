
import numpy as np
import sys
sys.path.append('../')

#import example
#if __name__ == '__main__':
#    xx = example.add(3, 4)
#    print(xx)

import SparseSolvPy

if __name__ == '__main__':
    temp_mat = [
        [10.0, 0.0, -0.1, 1.5, 0.01],
        [0.0, 15.6, 0.0, 0.0, 0.0],
        [-0.1, 0.0, 30.5, 0.0, 1.25],
        [1.5, 0.0, 0.0, 24.55, 0.0],
        [0.01, 0.0, 1.25, 0.0, 19.8]
    ]
    mat = SparseSolvPy.SparseMat(5)
    for i in range(5):
        for j in range(5):
            fval = abs(temp_mat[i][j])
            if fval > 1.0e-9:
                mat.add(i, j, temp_mat[i][j])
    mat.fix(False)
    #mat.print()
    mat.printMat("a.csv")

    mat2 = SparseSolvPy.SparseMatC(2)
    mat2.add(0, 0, 1.0)
    mat2.add(1, 5, 1.0j)
    mat2.fix(False)
    #mat2.print()
    mat2.printMat("a2.csv")

    gyo = [0, 0, 1, 1, 2, 3]
    retu =[0, 3, 0, 1, 5, 1]
    val =[0.1, 1.5, 2.5, -0.1, -5.0, 21.5]
    mat3 = SparseSolvPy.SparseMat(4, gyo, retu, val)
    mat3.printMat("a3.csv")


    # Set sparse solver
    solver = SparseSolvPy.MatSolvers()
    # save log
    solver.setSaveLog(True)

    #define vecB
    vec_b = np.asarray([1.0, 0.5, 0.0, 0.0, 0.0])
    results = np.asarray([0.0]*5)
    
    #solve !
    #bl = solver.solveICMRTR_py(5, 1.0e-6, 10000, 1.05, mat, vec_b, results, True)
    bl = solver.solveSGSMRTR_py(5, 1.0e-6, 10000, mat, vec_b, results, True)
    #bl = solver.solveICCGwithABMC_py(5, 1.0e-6, 10000, 1.05, mat, vec_b, results, 2, 4, True)
    print(results)

    # get residual log
    log = solver.getResidualLog_py()
    print(log)


#end

