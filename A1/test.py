import numpy as np
outerSide = 0.1
coreWidth = 0.04
coreHeight = 0.02
coreVoltage = 15.0
minres = 10**(-5)
h = 0.02
def createZeros(nRow, nCol):
    matrix = [[0 for i in range(nCol)] for j in range(nRow)]
    res = np.array(matrix, dtype=float)
    return res

def generateMesh(h):
    nRow = (int)(outerSide/h)+ 1
    nCol = (int)(outerSide/h) + 1
    mesh = createZeros(nRow,nCol)
    for i in range(nRow):
        for j in range(nCol):
            if(j<=(int(coreWidth/h))and i<= (int(coreHeight/h))):
                mesh[i][j] = coreVoltage
    return mesh


n = 5
for i in range(n):
    print(i)


