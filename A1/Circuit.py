from Matrix import matrix
#Author@Yi Zhu
#ID@260716006

m = matrix()

class Circuit(object):

    # Solve node voltages for a circuit read from txt file
    def solveCircuit(self, filename):
        file = open(filename, "r")
        circuit = file.readlines()
        J = list(map(float, circuit[0].split("\n")[0].split(" ")))
        R = list(map(float, circuit[1].split("\n")[0].split(" ")))
        E = list(map(float, circuit[2].split("\n")[0].split(" ")))
        A = []
        for line in circuit[4:]:
            A.append(list(map(float, line.split("\n")[0].split(" "))))
        Y = [[0 for x in range(len(R))] for y in range(len(R))]
        for i in range(len(Y)):
            Y[i][i] = 1 / R[i]
        At = m.matrixTranspose(A)
        AY = m.matrixMultiplication(A, Y)
        AYAT = m.matrixMultiplication(AY, At)
        YE = m.matrixVectorMultiplication(Y, E)
        JYE = m.vectorSubtraction(J, YE)
        b = m.matrixVectorMultiplication(A, JYE)
        res = m.solveMatrix(AYAT, b)
        return res

    # Find equivalent resistance for a N x 2N mesh circuit
    def findEqualRes(self, N, sR):
        A, J, R, E = self.generateNetwork(N, sR, 1)
        Y = [[0 for i in range(len(R))] for y in range(len(R))]
        for i in range(len(Y)):
            Y[i][i] = 1 / R[i]
        At = m.matrixTranspose(A)
        AY = m.matrixMultiplication(A, Y)
        AYAT = m.matrixMultiplication(AY, At)
        YE = m.matrixVectorMultiplication(Y, E)
        JYE = m.vectorSubtraction(J, YE)
        b = m.matrixVectorMultiplication(A, JYE)
        nodeVoltage = m.solveMatrix(AYAT, b)
        res = nodeVoltage[0] * sR / (1 - nodeVoltage[0])
        return res

    # Find equivalent resistance for a N x 2N mesh circuit use sparse matrix decomposition
    def sparseFindReq(self, N, sR):
        A, J, R, E = self.generateNetwork(N, sR, 1)
        Y = [[0 for i in range(len(R))] for y in range(len(R))]
        for i in range(len(Y)):
            Y[i][i] = 1 / R[i]
        band = N + 2
        At = m.matrixTranspose(A)
        AY = m.matrixMultiplication(A, Y)
        AYAT = m.matrixMultiplication(AY, At, )
        YE = m.matrixVectorMultiplication(Y, E)
        JYE = m.vectorSubtraction(J, YE)
        b = m.matrixVectorMultiplication(A, JYE)
        nodeVoltage = m.sparseSolveMatrix(AYAT, b, band)
        res = nodeVoltage[0] * sR / (1 - nodeVoltage[0])
        return res

    # Generate the X x 2N matrix with each resistor equals to 1k ohm
    def generateNetwork(self, N, sR, sE):
        node = (N + 1) * (2 * N + 1)
        branch = (N + 1) * 2 * N + (2 * N + 1) * N
        J = [0] * (branch + 1)
        R = [sR] * (branch + 1)
        E = [0] * (branch + 1)
        E[0] = sE
        A = [[0] * (branch + 1) for i in range(node - 1)]

        for i in range(node - 1):
            for j in range(branch + 1):
                row = i % (N + 1)
                column = i // (N + 1)
                if (j == 0 and i == 0):
                    A[i][j] = -1
                elif (j == (column * (2 * N + 1) + N + row + 1) and column < 2 * N):
                    A[i][j] = 1
                elif (j == (column * (2 * N + 1) + row + 1) and row < N):
                    A[i][j] = 1
                elif (j == ((column - 1) * (2 * N + 1) + N + row + 1) and column > 0):
                    A[i][j] = -1
                elif (j == (column * (2 * N + 1) + row) and row > 0):
                    A[i][j] = -1

        return A, J, R, E
