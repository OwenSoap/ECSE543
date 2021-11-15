import math, copy
#Author@Yi Zhu
#ID@260716006

class matrix(object):

    # Method to solve Ax=b using Choleski and fwd/bwd elimination
    def solveMatrix(self, matrix, b):
        self.checkSym(matrix)
        self.checkDet(matrix)
        L = self.choleskiDecompose(matrix)
        Y = self.forwardElm(L, b)
        Lt = self.matrixTranspose(L)
        X = self.backElm(Lt, Y)
        return X

    # Method to solve Ax=b using Choleski and fwd/bwd elimination and sparse matrix
    def sparseSolveMatrix(self, matrix, b, band):
        self.checkSym(matrix)
        self.checkDet(matrix)
        L = self.sparseCholeskiDecompose(matrix, band)
        Y = self.forwardElm(L, b)
        Lt = self.matrixTranspose(L)
        X = self.backElm(Lt, Y)
        return X

    # Choleski decomposition using look ahead method return L
    def choleskiDecompose(self, matrix):
        A = copy.deepcopy(matrix)
        n = len(A)
        L = [[0.0] * n for i in range(n)]
        for j in range(n):
            if A[j][j] <0:
                exit('Input matrix must be a positive-definite matrix!')
            L[j][j] = math.sqrt(A[j][j])
            for i in range(j + 1, n):
                L[i][j] = A[i][j] / L[j][j]
                for k in range(j + 1, i + 1):
                    A[i][k] = A[i][k] - L[i][j] * L[k][j]
        return L

    # Choleski decomposition using look ahead method and sparse matrix return L
    def sparseCholeskiDecompose(self, matrix, b):
        A = copy.deepcopy(matrix)
        n = len(A)
        L = [[0.0] * n for i in range(n)]
        for j in range(n):
            if A[j][j] < 0:
                exit('Input matrix must be a positive-definite matrix!')
            L[j][j] = math.sqrt(A[j][j])
            for i in range(j + 1, min(j + 1 + b, n)):
                L[i][j] = A[i][j] / L[j][j]
                for k in range(j + 1, min(j + 1 + b ,i + 1)):
                    A[i][k] = A[i][k] - L[i][j] * L[k][j]
        return L

    # Forward elimination return Y
    def forwardElm(self, lMatrix, bVector):
        L = copy.deepcopy(lMatrix)
        b = copy.deepcopy(bVector)
        n = len(b)
        Y = []
        for j in range(n):
            b[j] = b[j]/L[j][j]
            Y.append(b[j])
            for i in range(j+1, n):
                b[i] = b[i] - L[i][j]*b[j]
        return Y

    # Backward elimination return X
    def backElm(self, LtMatrix, yVector):
        Lt = copy.deepcopy(LtMatrix)
        Y = copy.deepcopy(yVector)
        n = len(Y)
        X = []
        for i in range(n)[::-1]:
            sum = 0
            for j in range(n)[:i:-1]:
                sum += X[n - j - 1] * Lt[i][j]
            X.append((Y[i] - sum) / Lt[i][i])
        return X[::-1]

    # Check if the matrix is symmetric
    def checkSym(self, matrix):
        n = len(matrix)
        for i in range(n):
            for j in range(i + 1, n):
                if matrix[i][j] != matrix[j][i]:
                    exit('Input matrix must be a symmetric matrix!')

    # Check if the determinant of the matrix is zero
    def checkDet(self, matrix):
        n = len(matrix)
        det = 1
        for i in range(n):
            det *= matrix[i][i]
        if det <= 0:
            exit('Input matrix must be a positive definite matrix!')

    # Check if A * x equals b
    def checkSol(self, A, X, b):
        n = len(A)
        for i in range(n):
            res = self.matrixVectorMultiplication(A,X)
            if abs(b[i] - res[i]) > 0.01:
                return False
        return True

    # Matrix transpose return the transpose of a matrix
    def matrixTranspose(self, matrix):
        A = copy.deepcopy(matrix)
        nRow = len(A)
        nCol = len(A[0])
        res = []
        for i in range(nCol):
            row = []
            for j in range(nRow):
                row.append(A[j][i])
            res.append(row)
        return res

    # Matrix multiply by a vector return the result vector
    def matrixVectorMultiplication(self, A, b):
        res = []
        nRow = len(A)
        nCol = len(b)
        for i in range(nRow):
            sum = 0
            for j in range(nCol):
                sum += A[i][j] * b[j]
            res.append(sum)
        return res

    # Matrix multiply by another matrix return the result matrix
    def matrixMultiplication(self, A, B):
        nRow = len(A)
        nCol = len(B[0])
        nB = len(B)
        res = [0.0] * nRow
        for i in range(nRow):
            res[i] = [0] * nCol
        for i in range(nRow):
            for j in range(nCol):
                for k in range(nB):
                    res[i][j] += A[i][k] * B[k][j]
        return res

    # Vector multiply by another matrix return the result vector
    def vectorMatrixMultiplication(self, b, A):
        res = []
        nRow = len(A)
        nCol = len(b)
        for i in range(nRow):
            sum = 0
            for j in range(nCol):
                sum += b[j] * A[i][j]
            res.append(sum)
        return res

    # Two vectors adding up return the result vector
    def vectorAddition(self, a, b):
        n = len(a)
        res = []
        for i in range(n):
            res.append(a[i] + b[i])
        return res

    # Two vectors subtracting return result vector
    def vectorSubtraction(self, a, b):
        n = len(a)
        res = []
        for i in range(n):
            res.append(a[i] - b[i])
        return res

    # Two vectors multiplying return result value
    def vectorMultiplication(self, a, b):
        nb = len(b)
        res = 0
        for i in range(nb):
            res += a[i] * b[i]
        return res

    # Two matrix subtraction return result matrix
    def matrixSubstract(self, A, B):
        res = []
        nRow = len(A)
        nCol = len(A[0])
        for i in range(nRow):
            row = []
            for j in range(nCol):
                row.append(A[i][j]-B[i][j])
            res.append(row)
        return res

    # 2x2 matrix inverse return the inverse of a 2x2 matrix A
    def matrixInverse(self, A):
        det = 1.0 / (A[0][0] * A[1][1] - A[1][0] * A[0][1])
        temp = A[0][0]
        A[0][0] = A[1][1] * det
        A[1][1] = temp * det
        A[1][0] *= -1.0 * det
        A[0][1] *= -1.0 * det

        return A



