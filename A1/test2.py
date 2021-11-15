import time

from Matrix import matrix
from Circuit import Circuit
m = matrix()
c = Circuit()

#TEST CHEOLSKI
def cholskiTest():
    A = [[3,1,1],[1,3,2],[1,2,5]]
    r = m.choleskiDecompose(A)
    rt = m.matrixTranspose(r)
    ans = m.matrixMultiplication(r, rt)
    for i in range(len(ans)):
        for j in range(len(ans)):
            ans[i][j] = round(ans[i][j],2)
    print(str(r) + "\n")
    print(str(rt) + "\n")
    print(str(ans) + "\n")

#TEST SOLVE MATRIX
def solveMatrixTest():
    A = [[3,1,1],[1,3,2],[1,2,5]]
    b = [8,13,20]
    x = [1,2,3]

    r = m.solveMatrix(A,b)
    for i in range(len(r)):
        r[i] = round(r[i],2)
    ans = m.checkSol(A,r,b)
    print("Matrix A is "+ str(A) + "\n" + "Vector b is "+ str(b) + "\n" + "Expect x value to be " + str(x) + "\n" +
            "Calculated result x is " + str(r) + "\n" + "The calculation is " + str(ans))

#TEST SOLVE CIRCUITS NODE VOLTAGE
def solveCircuitsVoltage():
    res= c.solveCircuit("testCircuit4.txt")
    print("Node voltage found for testCircuit1: " + str(res))

def generateCiurcuit():
    N = 1
    sR = 1000
    A, J, R, E = c.generateNetwork(N, sR, 1)
    print("J: " + str(J))
    print("R: " + str(R))
    print("E: " + str(E))
    print("A: " + str(A))

#TEST FIND EQUIVLENT RESISTOR
def findEqualRes():
    N = 4
    nodeNum = (N + 1) * (2 * N + 1)
    branchNum = (N + 1) * 2 * N + (2 * N + 1) * N
    start = time.time()
    ans = c.findEqualRes(N,1)
    runTime = (time.time() - start)
    res = ans * 1000
    #print("N = " + str(N))
    #print("Node number: " + str(nodeNum))
    #print("Mesh number: " + str(branchNum))
    print("Equivalent resistor is " + str(res) + "ohm" + "\n"+"Runtime is: ")
    print(str(time.time() - start)+"s")

def sparseFindEqualRes():
    N = 4
    nodeNum = (N + 1) * (2 * N + 1)
    branchNum = (N + 1) * 2 * N + (2 * N + 1) * N
    start = time.time()
    ans = c.sparseFindReq(N, 1)
    runTime = (time.time() - start)
    res = ans * 1000
    # print("N = " + str(N))
    # print("Node number: " + str(nodeNum))
    # print("Mesh number: " + str(branchNum))
    print("Equivalent resistor is " + str(res) + "ohm" + "\n"+"Runtime is: ")
    print(str(runTime)+"s")


class test(object):

    if __name__=="__main__":
        solveCircuitsVoltage()