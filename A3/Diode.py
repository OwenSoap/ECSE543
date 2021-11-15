import math
from Matrix import matrix

E = 0.22
R = 500
Isa = 0.6e-6
Isb = 1.2e-6
ktq = 25e-3

m = matrix()

# Author@Yi Zhu
# ID@260716006

class Diode(object):

    # Function solves nonlinear equation using Newton-Raphson iteration method
    def newtonRaphson(self, vn, maxerror):
        count = 0
        v1 = vn[0]
        v2 = vn[1]

        f1 = (E - v1) / R - Isa * (math.exp((v1 - v2) / ktq) - 1.0)
        f2 = (E - v1) / R - Isb * (math.exp(v2 / ktq) - 1.0)
        f = [f1, f2]

        J = [[0 for x in range(2)] for y in range(2)]
        J[0][0] = (-1 / R) - (Isa / ktq) * (math.exp((v1 - v2) / ktq))
        J[0][1] = (Isa / ktq) * (math.exp((v1 - v2) / ktq))
        J[1][0] = (-1 / R)
        J[1][1] = -1 * (Isb / ktq) * (math.exp(v2 / ktq))

        Jinv = m.matrixInverse(J)
        Jinvf = m.matrixVectorMultiplication(Jinv, f)
        vn = m.vectorSubtraction(vn, Jinvf)
        error = [abs(z) for z in Jinvf]
        print("Iteration: " + str(count) + " v1: " + str(vn[0]) + "v v2: " + str(vn[1]) + "v f: " + str(f))

        while (abs(max(error)) > maxerror):
            count += 1
            v1 = vn[0]
            v2 = vn[1]

            f1 = (E - v1) / R - Isa * (math.exp((v1 - v2) / ktq) - 1.0)
            f2 = (E - v1) / R - Isb * (math.exp(v2 / ktq) - 1.0)
            f = [f1, f2]

            J = [[0 for x in range(2)] for y in range(2)]
            J[0][0] = (-1 / R) - (Isa / ktq) * (math.exp((v1 - v2) / ktq))
            J[0][1] = (Isa / ktq) * (math.exp((v1 - v2) / ktq))
            J[1][0] = (-1 / R)
            J[1][1] = -1 * (Isb / ktq) * (math.exp(v2 / ktq))

            Jinv = m.matrixInverse(J)
            Jinvf = m.matrixVectorMultiplication(Jinv, f)
            vn = m.vectorSubtraction(vn, Jinvf)
            error = [abs(z) for z in Jinvf]
            print("Iteration: " + str(count) + " v1: " + str(vn[0]) + "v v2: " + str(vn[1]) + "v f: " + str(f))

        return f, Jinvf, vn, count
