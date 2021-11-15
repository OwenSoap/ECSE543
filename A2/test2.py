from Capacitance import Capacitance
from Matrix import matrix
from ConjugateGradient import ConjugateGradient

c = Capacitance()
m = matrix()
cg = ConjugateGradient()
U = c.readResult("result.dat")
s = c.S
ufix = [15, 15, 15, 15]
enery = c.generateMeshEnergy(U, s, ufix)
#res = c.capacitance(enery)
#print("Energy is " + str(enery) + " J")
#print("Capacitance is " + str(res) + " F")

A = cg.generateMeshMatrix()
b = cg.generateBVector()
At = m.matrixTranspose(A)
AtA = m.matrixMultiplication(At, A)
Atb = m.matrixVectorMultiplication(At,b)
resChe = m.solveMatrix(AtA,Atb)
res = cg.ConjugateGradient(AtA,Atb)
#print(A)
#print(AtA)

#for i in range(len(resChe)):
#   resChe[i]= round(resChe[i], 4)
#for i in range(len(resCon)):
#   resCon[i] = round(resCon[i], 4)
#print(AtA)
#print(Atb)
#print("The result vector x solved by using the Choleski Decomposition is: \n " + str(resChe))
res1 = res[1]
res2 = res[2]
for i in range(len(res1)):
  res1[i] = round(res1[i], 4)

for i in range(len(res2)):
  res2[i] = round(res2[i], 4)
print(res1)
print(res2)
#print("The result vector x solved by using the Conjugate Gradient method is: \n " + str(resCon))
#print("The potential at (x,y) = (0.06, 0.04) solved by using the Choleski Decomposition is: \n "+str(resChe[11])+"v")
#print("The potential at (x,y) = (0.06, 0.04) solved by using the Conjugate Gradient method is: \n "+str(res1[11])+"v")
