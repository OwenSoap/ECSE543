from Matrix import matrix
from math import sqrt
import copy

# Author@Yi Zhu
# ID@260716006

m = matrix()


class ConjugateGradient(object):
    def __init__(self):
        pass

    # Generate a mesh Matrix A containing voltages at free node using five-points method
    def generateMeshMatrix(self):
        A = [[0.0] * 19 for i in range(19)]
        for i in range(19):
            A[i][i] = -4
            if (i < 10):
                if (i < 8):
                    A[i][i + 5] = 1
                if (i != 0 and i != 5):
                    A[i][i - 1] = 1
                if (i + 1 == 5 or i + 1 == 10):
                    A[i][i - 1] += 1
                else:
                    A[i][i + 1] = 1
                if i > 4:
                    A[i][i - 5] = 1

            elif (i < 13):
                A[i][i - 5] = 1
                A[i][i + 3] = 1
                if (i != 10):
                    A[i][i - 1] = 1
                if (i != 12):
                    A[i][i + 1] = 1
            else:
                A[i][i - 3] = 1
                if (i < 16):
                    A[i][i + 3] = 1
                else:
                    A[i][i - 3] += 1
                if (i != 15) and (i != 18):
                    A[i][i + 1] = 1
                if (i != 13) and (i != 16):
                    A[i][i - 1] = 1
        return A

    # Generate a vector b containing voltages at fixed node
    def generateBVector(self):
        b = [0.0] * 19
        for i in range(19):
            if i == 8 or i == 9 or i == 12 or i == 15 or i == 18:
                b[i] = -15
        return b

    # Function solve Ax=b using Conjugate Gradient method and calculate infinity norm and 2-norm
    # Return result x, infinity norm and 2-norm
    def ConjugateGradient(self, A, b):
        infNorm = []
        twoNorm = []
        nCol = len(b)
        x = [0 for i in range(nCol)]
        r = m.vectorSubtraction(b, m.matrixVectorMultiplication(A, x))
        p = copy.deepcopy(r)

        infiniteNorm = max([abs(j) for j in r])
        infNorm.append(infiniteNorm)
        res = 0
        norm = 0
        for k in r:
            res += k ** 2
            norm = sqrt(res)
        twoNorm.append(norm)

        for i in range(nCol):
            alpha = (m.vectorMultiplication(p, r)) / (m.vectorMultiplication(m.vectorMatrixMultiplication(p, A), p))
            x = m.vectorAddition(x, [alpha * j for j in p])
            r = m.vectorSubtraction(b, m.matrixVectorMultiplication(A, x))
            beta = -1 * ((m.vectorMultiplication(p, m.matrixVectorMultiplication(A, r))) / (
                m.vectorMultiplication(p, m.matrixVectorMultiplication(A, p))))
            p = m.vectorAddition(r, [beta * j for j in p])

            infiniteNorm = max([abs(k) for k in r])
            infNorm.append(infiniteNorm)
            res = 0
            norm = 0
            for l in r:
                res += l ** 2
                norm = sqrt(res)
            twoNorm.append(norm)

        return x, infNorm, twoNorm
