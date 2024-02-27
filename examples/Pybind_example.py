
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import SparseSolvPy
import numpy as np

if __name__ == '__main__':
    # Sample for simple use
    mat = SparseSolvPy.SparseMat(10)
    mat.add(0, 0, 1.0)
    mat.fix(False)
    #mat.print()
    mat.printMat("a.csv")
    # Sample for simple use (Complex)
    mat2 = SparseSolvPy.SparseMatC(2)
    mat2.add(0, 0, 1.0)
    mat2.add(1, 5, 1.0j)
    mat2.fix(False)
    #mat2.print()
    mat2.printMat("a2.csv")

    #Sample (direct init. by sparse data)
    gyo = [0, 0, 1, 1, 2, 3]
    retu =[0, 3, 0, 1, 5, 1]
    val = np.asarray([0.1, 1.5, 2.5, -0.1, -5.0, 21.5])
    mat3 = SparseSolvPy.SparseMat(6, gyo, retu, val)
    mat3.printMat("a3.csv")


    #sample for ICCG
    mat_test = [
        [ 1.5,  -0.5, 0.0,  0.0,   1.6],
        [-0.5,   5.2, 0.0, -0.1,   2.5],
        [ 0.0,   0.0, 8.5,  0.0,   0.0],
        [ 0.0,  -0.1, 0.0,  7.7,   0.6],
        [ 1.6,   2.5, 0.0,  0.6,  10.6]
    ]
    mat4 = SparseSolvPy.SparseMat(5)
    for i in range(5):
        for j in range(5):
            val = abs(mat_test[i][j])
            if val > 1.0e-12:
                mat4.add(i, j, mat_test[i][j])
            #end
        #end
    #end
    mat4.fix(False)
    mat4.printMat("a4.csv")

    vecB = np.asarray([1.0, 1.0, 5.0, 0.1, 5.0])
    results = np.asarray([0.0]*5)

    log = SparseSolvPy.MatSolversICCG.solveICCGPy(5, 1.0e-6, 10000, 1.05, mat4, vecB, results, True)
    print(results)
    for i in range(5):
        temp = 0.0        
        for j in range(5):
            temp += mat_test[i][j] * results[j]
        #end
        print(i, temp, vecB[i])
    #end
    print("log is ")
    for arg in log:
        print(arg)
#end
    