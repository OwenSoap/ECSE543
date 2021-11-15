from Matrix import matrix

# Author@Yi Zhu
# ID@260716006

matrix = matrix()


class Capacitance(object):
    def __init__(self):
        pass

    S = [[1, -0.5, 0, -0.5], [-0.5, 1, -0.5, 0], [0, -0.5, 1, -0.5], [-0.5, 0, -0.5, 1]]

    # Function generating the energy matrix using voltage at free nodes, fixed nodes and Sglobal
    def generateMeshEnergy(self, U, Sg, uFix):
        n = 0
        m = 0
        sum = 0
        while (n < 28):
            if (n == 5 or n == 11 or n == 17 or n == 18 or n == 23):
                n += 1
            if (n < 17):
                uV = self.part1(U, n)
                w = matrix.matrixVectorMultiplication(Sg, uV)
                e = matrix.vectorMultiplication(uV, w)
                sum += 0.5 * e
                n += 1
                m += 1
            elif (n > 18 and n < 28):
                uV = self.part2(U, n)
                w = matrix.matrixVectorMultiplication(Sg, uV)
                e = matrix.vectorMultiplication(uV, w)
                sum += 0.5 * e
                n += 1
                m += 1

        wFix = matrix.matrixVectorMultiplication(Sg, uFix)
        eFix = matrix.vectorMultiplication(uFix, wFix)
        sum += 2 * 0.5 * eFix
        sum = sum * 8.85 * pow(10, -12)
        return sum

    # Stamp free node voltage from U to nodes number <17
    def part1(self, U, n):
        uVector = []
        a = U[n]
        b = U[n + 1]
        c = U[n + 7]
        d = U[n + 6]
        uVector.append(a)
        uVector.append(b)
        uVector.append(c)
        uVector.append(d)
        return uVector

    # Stamp free node voltage from U to nodes number <28
    def part2(self, U, n):
        uVector = []
        a = U[n]
        b = U[n + 1]
        c = U[n + 6]
        d = U[n + 5]
        uVector.append(a)
        uVector.append(b)
        uVector.append(c)
        uVector.append(d)
        return uVector

    # Calculate capacitance per unit length using formula C=2E/V^2
    def capacitance(self, energy):
        return 2 * 4 * energy / (15 * 15)

    # Read the potential from output file of SIMPLE_2D method
    def readResult(self, filename):
        U = []
        file = open(filename, "r")
        potential = file.readlines()
        for row in potential:
            U.append(float(row.split("\n")[0].split(" ")[-1]))
        return U
