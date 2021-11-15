from Matrix import matrix
m = matrix()
"""
from sympy import symbols, expand, diff, lambdify
def full_lagrange_interpolate(X, Y):
    # lagrange dimension decided by length of X and Y
    n = len(X)
    x = symbols('x')
    Ls = []
    # find Fj for each j in n
    for j in range(n):
        X_sub = X[:j]+X[j+1:]
        F = 1
        for xr in X_sub:
            F *= x-xr
        F_f = lambdify(x, F)
        L = expand(F/F_f(X[j]))
        Ls.append(L)
    y = 0
    for j, L in enumerate(Ls):
        y += Y[j]*L
    return expand(y)

B = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4]
res = full_lagrange_interpolate(B,H)
print(res)
"""
J = [[0 for x in range(2)] for y in range(2)]
J[0][0] = 5
J[0][1] = 2
J[1][0] = -7
J[1][1] = -3
print(J)
iJ=m.matrixInverse(J)
print(iJ)
