import numpy as np

outerSide = 0.1
coreWidth = 0.04
coreHeight = 0.02
coreVoltage = 15.0
minres = 10**(-5)

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

def SOR(mesh,h,w):
    for i in range (len(mesh) - 1):
        for j in range(len(mesh[i]) - 1):
            if (j >=int(coreWidth/h) or i >= int(coreHeight/h)):
                a = mesh[i - 1][j]
                b = mesh[i][j - 1]
                if((i-2)<=0):
                    a = mesh[i][j]
                if((j-1<0)):
                    b = mesh[i][j]
                mesh[i][j] = (1-w)*mesh[i][j]+(w/4)*(a+b+mesh[i][j+1]+mesh[i+1][j])
    return mesh

def maxResidual(mesh,h,minRes):
    nodeI = (int)(outerSide/h)
    nodeJ = (int)(outerSide/h)
    max = 0
    for i in range(0,nodeI):
        for j in range(0,nodeJ):
            if (j >= int(coreWidth / h) or i >= int(coreHeight / h)):
                if (j >= int(coreWidth / h) or i >= int(coreHeight / h)):
                    a = mesh[i - 1][j]
                    b = mesh[i][j - 1]
                    if ((i - 1) < 0):
                        a = mesh[i][j]
                    if ((j - 1 < 0)):
                        b = mesh[i][j]
                residual = a+b+mesh[i][j+1]+mesh[i+1][j] - 4*mesh[i][j]
                residual = abs(residual)
                if(residual>max):
                    max = residual
    if(max>=minRes):
        return True
    else:
        return False

def SORIter(mesh,h,w,minRes):
    iteration = 0
    while(maxResidual(mesh,h,minRes)):
        mesh = SOR(mesh,h,w)
        iteration+=1
    result = {
        'mesh': mesh,
        'iteration': iteration
    }
    return result

def Jacobian(mesh,h):
    nodeI = (int)(outerSide / h)
    nodeJ = (int)(outerSide / h)
    for i in range (0,nodeI):
        for j in range(0,nodeJ):
            if (i>=(int)(coreWidth/h) or j >= (int)(coreHeight/h)):
                a = mesh[i - 1][j]
                b = mesh[i][j - 1]
                if ((i - 1) < 0):
                    a = mesh[i][j]
                if ((j - 1 < 0)):
                    b = mesh[i][j]
                mesh[i][j]=(a+b+mesh[i][j+1]+mesh[i+1][j])/4
    return mesh

def jacIter(mesh,h,minRes):
    iteration = 0
    while(maxResidual(mesh,h,minRes)):
        mesh = Jacobian(mesh,h)
        iteration+=1
    result = {
        'mesh':mesh,
        'iteration':iteration
    }
    return result


def getVoltage(mesh,x,y,h):
    xNode = (int)(x/h)
    yNode = (int)(y/h)
    return mesh[xNode-1][yNode-1]

"""
#b
# for i in range(0,10):
#     w = 1 + 0.1*i+0.02
#     h = 0.02
#     mesh = generateMesh(h)
#     result = SORIter(mesh,h,w,minres)
#     finMesh = result['mesh']
#     voltage = getVoltage(finMesh,0.06,0.04,h)
#     finalIter = result['iteration']
#     print('for w = '+str(w)+' the voltage is '+str(voltage)+' iteration is '+str(finalIter))

h = 0.02
for i in range(1,11):
    h*=0.8
    w = 1.5
    mesh = generateMesh(h)
    result = SORIter(mesh,h,w,minres)
    finMesh = result['mesh']
    voltage = getVoltage(finMesh,0.06,0.04,h)
    finalIter = result['iteration']
    print('for h = '+str(h)+' the voltage is '+str(voltage)+' iteration is '+str(finalIter))
d
h = 0.02
for i in range(1,11):
    h*=0.8
    w = 1.5
    mesh = generateMesh(h)
    result = jacIter(mesh,h,minres)
    finMesh = result['mesh']
    voltage = getVoltage(finMesh,0.06,0.04,h)
    finalIter = result['iteration']
    print('for h = '+str(h)+' the voltage is '+str(voltage)+' iteration is '+str(finalIter))
"""
h = 0.02
w = 0.00001
mesh = generateMesh(h)
res = SOR(mesh,h,w)
#for i in range(len(res[0])):
    #for j in range(len(res)):
        #res[i][j] = round(res[i][j],10)
print(res)