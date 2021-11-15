import numpy as np
import copy

# Author@Yi Zhu
# ID@260716006

outerWidth = 0.1
innerWidth = 0.04
innerHeight = 0.02
innerVoltage = 15.0
minResidual = 10 ** (-5)


# Create a zero matrix in matrix form
def createZeros(nRow, nCol):
    matrix = [[0 for i in range(nCol)] for j in range(nRow)]
    res = np.array(matrix, dtype=float)
    return res


# Generate mesh matrix with node voltage
def generateMesh(h):
    nRow = int(outerWidth / h) + 1
    nCol = int(outerWidth / h) + 1
    mesh = createZeros(nRow, nCol)
    for i in range(nRow):
        for j in range(nCol):
            if (j <= (int(innerWidth / h)) and i <= (int(innerHeight / h))):
                mesh[i][j] = innerVoltage
    return mesh


# Solve mesh matrix using SOR method, return new mesh and iteration number
def SOR(mesh, h, w, residual):
    nRow = int(outerWidth / h)
    nCol = int(outerWidth / h)
    x = int(innerWidth / h) + 1
    y = int(innerHeight / h) + 1
    it = 0
    while (underResidual(mesh, h, residual)):
        for i in range(nRow):
            for j in range(nCol):
                if (i >= y or j >= x):
                    a = mesh[i - 1][j]
                    b = mesh[i][j - 1]
                    if (i == 0):
                        a = mesh[i][j]
                    if (j == 0):
                        b = mesh[i][j]
                    mesh[i][j] = (1 - w) * mesh[i][j] + (w / 4) * (a + b + mesh[i][j + 1] + mesh[i + 1][j])
        it += 1
    return mesh, it


# Solve mesh matrix using Jacobi method, return new mesh and iteration number
def JacoBi(mesh, h, residual):
    nRow = int(outerWidth / h)
    nCol = int(outerWidth / h)
    x = int(innerWidth / h) + 1
    y = int(innerHeight / h) + 1
    it = 0
    while (underResidual(mesh, h, residual)):
        m = copy.deepcopy(mesh)
        for i in range(nRow):
            for j in range(nCol):
                if (i >= y or j >= x):
                    a = m[i - 1][j]
                    b = m[i][j - 1]
                    if (i == 0):
                        a = mesh[i][j]
                    if (j == 0):
                        b = mesh[i][j]
                    mesh[i][j] = (1 / 4) * (a + b + mesh[i][j + 1] + mesh[i + 1][j])
        it += 1
    return mesh, it


# Solve mesh matrix using SOR method, return new mesh and iteration number
def nonUniform(mesh, h, w, minResidual, x_space, y_space):
    nRow = int(outerWidth / h)
    nCol = int(outerWidth / h)
    x = int(innerWidth / h) + 1
    y = int(innerHeight / h) + 1
    it = 0
    max = 0
    inResidual = True
    while (inResidual):
        for i in range(nRow):
            for j in range(nCol):
                if (i >= y or j >= x):
                    a1 = abs(x_space[i] - x_space[i - 1])
                    a2 = abs(x_space[i + 1] - x_space[i])
                    b1 = abs(y_space[j] - y_space[j - 1])
                    b2 = abs(y_space[j + 1] - y_space[j])
                    sum = (mesh[i - 1][j] / (a1 * (a1 + a2)) + mesh[i + 1][j] / (a2 * (a1 + a2)) + mesh[i][j - 1] / (b1 * (b1 + b2)) + mesh[i][j + 1] / (b2 * (b1 + b2))) / (1 / (a1 * (a1 + a2)) + 1 / (a2 * (a1 + a2)) + 1 / (b1 * (b1 + b2)) + 1 / (b2 * (b1 + b2)))
                    vol = (1 - w) * mesh[i][j] + w * sum
                    residual = vol - mesh[i][j]
                    residual = abs(residual)
                    if (residual > max):
                        max = residual
                    mesh[i][j] = vol
        if (max < minResidual):
            inResidual = False

        it += 1
    return mesh, it


# Check whether residual is less than min residual, return boolean value true or false
def underResidual(mesh, h, minResidual):
    nRow = int(outerWidth / h)
    nCol = int(outerWidth / h)
    x = int(innerWidth / h) + 1
    y = int(innerHeight / h) + 1
    max = 0
    for i in range(nRow):
        for j in range(nCol):
            if (i >= y or j >= x):
                a = mesh[i - 1][j]
                b = mesh[i][j - 1]
                if (i == 0):
                    a = mesh[i][j]
                if (j == 0):
                    b = mesh[i][j]
                residual = a + b + mesh[i][j + 1] + mesh[i + 1][j] - 4 * mesh[i][j]
                residual = abs(residual)
                if (residual > max):
                    max = residual
    if (max >= minResidual):
        return True
    else:
        return False


# Get node voltage of specific point
def getVoltage(mesh, x, y, h):
    row = int(x / h)
    col = int(y / h)
    return mesh[row][col]


x_space = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
y_space = [0, 0.02, 0.04, 0.05, 0.055, 0.06, 0.065, 0.07, 0.08, 0.085, 0.1]

h = 0.01
w = 1.4
mesh = generateMesh(h)
res, i = SOR(mesh, h, w, minResidual)
#res, i = SOR(mesh, h, w, minResidual)
#res, i = nonUniform(mesh, h, w, minResidual, x_space, y_space)
voltage = getVoltage(res, 0.06, 0.04, h)
#print("h: " + str(round(h, 4)) + " w: " + str(w))
print("Iteration: " + str(i) + " Node(0.06,0.04) voltage is: " + str(voltage) + "v")
