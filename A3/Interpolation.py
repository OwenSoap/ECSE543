from sympy import symbols, expand, lambdify, diff

# Author@Yi Zhu
# ID@260716006

class Interpolation(object):

    # Function operates interpolation using full-domain Lagrange polynomials
    def lagrange(sefl, X, Y):
        x = symbols('x')
        res = 0
        n = len(X)

        def l(i, n):
            lx = 1
            for j in range(n):
                if i == j:
                    continue
                xj = X[j]
                lx *= (x - xj) / (xi - xj)
            return lx

        for i in range(n):
            xi = X[i]
            yi = Y[i]
            res += yi * l(i, n)

        return expand(res)

    # Function operates interpolation using Cubic Hermite polynomials
    def cubicHermite(self, X, Y):
        x = symbols('x')
        n = len(X)
        res = 0
        Uj = []
        Vj = []

        def l(i, n):
            lx = 1
            for j in range(n):
                if i == j:
                    continue
                xj = X[j]
                lx *= (x - xj) / (xi - xj)
            return lx

        for i in range(n):
            xi = X[i]
            lx = l(i, n)
            ld = lambdify(x, diff(lx))
            U = (1 - 2 * ld(X[i]) * (x - X[i])) * (lx ** 2)
            V = (x - X[i]) * (lx ** 2)
            Uj.append(U)
            Vj.append(V)

        Yprime = []
        for i in range(n - 1):
            Yprime.append((Y[i + 1] - Y[i]) / (X[i + 1] - X[i]))

        Yprime.append(Y[-1] / X[-1])

        for j in range(n):
            res += Y[j] * Uj[j] + Yprime[j] * Vj[j]

        return expand(res)
