from Capacitance import Capacitance
from Matrix import matrix
from ConjugateGradient import ConjugateGradient
import copy, math
cg = ConjugateGradient()
c = Capacitance()
m = matrix()
A = cg.generateMeshMatrix()
b = cg.generateBVector()
At = m.matrixTranspose(A)
AtA = m.matrixMultiplication(At, A)
Atb = m.matrixVectorMultiplication(At,b)
res = m.solveMatrix(AtA,Atb)

# pTr/pTAp
def alpha(self, A, r, p):
    pTr = m.vectorMultiplication(p, r)
    pTA = m.vectorMatrixMultiplication(p, A)
    pTAp = m.vectorMultiplication(pTA, p)
    return float(pTr) / float(pTAp)

def residue(self, A, x, b):
    Ax = m.matrixVectorMultiplication(A, x)
    return m.vectorSubtraction(b, Ax)

def newGuess(self, x, alpha, p):
    alphap = [alpha * i for i in p]
    return m.vectorAddition(x, alphap)

# -pTAr/pTAp
def beta(self, A, r, p):
    Ar = m.matrixVectorMultiplication(A, r)
    pTAr = m.vectorMultiplication(p, Ar)
    Ap = m.matrixVectorMultiplication(A, p)
    pTAp = m.vectorMultiplication(p, Ap)
    return -1 * float(pTAr) / float(pTAp)

# newR + bp
def newOrientation(self, r, beta, p):
    bp = m.scaleVector(beta, p)
    return m.vectorAddition(r, bp)

def ConjugateGradient(A, b):
    nCol = len(b)
    x = [0 for i in range(nCol)]
    r = m.vectorSubtraction(b, m.matrixVectorMultiplication(A, x))
    p = copy.deepcopy(r)
    for i in range(nCol):
        alpha = (m.vectorMultiplication(p, r))/(m.vectorMultiplication(m.vectorMatrixMultiplication(p, A), p))
        x = m.vectorAddition(x, [alpha * i for i in p])
        r = m.vectorSubtraction(b, m.matrixVectorMultiplication(A, x))
        beta = -1 * ((m.vectorMultiplication(p, m.matrixVectorMultiplication(A, r)))/(m.vectorMultiplication(p, m.matrixVectorMultiplication(A, p))))
        p = m.vectorAddition(r, [beta * i for i in p])

    return x

A = ConjugateGradient(AtA,Atb)
print(A)
