from sympy import symbols, expand, diff, lambdify
import IronCore, Diode, Interpolation, math, Integral
ic = IronCore.IronCore()
dd = Diode.Diode()
ip = Interpolation.Interpolation()
i =  Integral.Integral()
""""
def Lagrange(X, Y):
    x = symbols('x')
    p = 0
    n = len(X)
    for i in range(n):
        xi = X[i]
        yi = Y[i]

        def l(i, n):
            lx = 1
            for j in range(n):
                if i == j:
                    continue
                xj = X[j]
                lx *= (x - xj) / float(xi - xj)
            return lx

        p += yi * l(i, n)
    return expand(p)

def Cubic_Hermite(X, Y):
    x = symbols('x')
    n = len(X)
    p = 0
    Us = []
    Vs = []
    for i in range(n):
        xi = X[i]

        def l(i, n):
            lx = 1
            for j in range(n):
                if i == j:
                    continue
                xj = X[j]
                lx *= (x - xj) / float(xi - xj)
            return lx

        lx = l(i, n)
        ld = lambdify(x, diff(lx))
        U = (1 - 2 * ld(X[i]) * (x - X[i]))*(lx**2)
        V = (x - X[i])*(lx**2)
        Us.append(U)
        Vs.append(V)

    b = []
    for j in range(n-1):
        b.append((Y[j+1] - Y[j]) / (X[j+1] - X[j]))
    b.append(Y[-1]/X[-1])

    for k in range(n):
        p += Y[k] * Us[k] + b[k] * Vs[k]

    return expand(p)

B = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4]

B1 = [0.0, 1.3, 1.4, 1.7, 1.8, 1.9]
H1 = [0.0, 540.6, 1062.8, 8687.4, 13924.3, 22650.2]
P = ip.lagrange(B, H)
P1 = ip.lagrange(B1,H1)
P2 = ip.cubicHermite(B1, H1)
print("Full-domain Lagrange Interpolation: ")
print(P)
print("Full-domain Lagrange Interpolation with certain B: ")
print(P1)
print("Cubic Hermite Polynomials Interpolation: ")
print(P2)

result = ic.newtonRaphson(0.0, 1e-6)
flux = result[0]
count = result[1]

print("Flux in the M19 iron core calculated by NR method is: ", round(flux, 10), " Wb")
print("Iteration: ", count)

result1 = ic.successiveSubstitution(0.0, 1e-6)
flux1 = result1[0]
count1 = result1[1]
print("Flux in the M19 iron core calculated by SS method is: ", round(flux1, 10), " Wb")
print("Iteration: ", count1)


v = [0, 0]

result = dd.newtonRaphson(v, 1e-6)
f = result[0]
j = result[1]
v = result[2]
c = result[3]
#print(f)
#print("Iteration: " + str(c) + " v1: " + str(v[0]) + " v2: " + str(v[1]) + " f: " + str(f))
# print(c)

"""
ans1 = math.sin(1)

for n in range(1, 21):
    res = i.gaussLegendreUni(math.cos, n, 0, 1)
    error = abs(ans1 - res)
    print("N = ", n, "  Result = ", res, "  Absolute Error = ", error)


ans2 = -1
for m in [(10 * d) for d in range(1, 21)]:
    res = i.gaussLegendreUni(math.log, m, 0, 1)
    error = abs(ans2 - res)
    print("N = ", m, "  Result = ", res, "  Absolute Error = ", error)


coord = [0.0, 0.03, 0.04, 0.08, 0.13, 0.22, 0.45, 0.68, 0.73, 0.91, 1]
ans2 = -1
res = i.gaussLegendreNonUni(math.log, coord)
error = abs(ans2 - res)
print( "Result = " + str(res) + "  Absolute Error = " + str(error))

