B = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
H = [0.0, 14.7, 36.5, 71.7, 121.4, 197.4, 256.2, 348.7, 540.6, 1062.8, 2318.0, 4781.8, 8687.4, 13924.3, 22650.2]
A = 0.0001
Lc = 0.3
La = 0.5
Ra = 3.9788935e7
NI = 8000

# Author@Yi Zhu
# ID@260716006

class IronCore(object):

    # Function solves nonlinear equation using Newton-Raphson iteration method
    def newtonRaphson(self, iguess, maxerror):
        count = 0
        psi = iguess
        H = self.H(psi)
        fi = Lc * H - NI
        while (abs(self.f(psi) / fi) > maxerror):
            psi -= self.f(psi) / (Ra + Lc * self.Hprime(psi) / A)
            count += 1
        return psi, count

    # Nonlinear function f
    def f(self, psi):
        f = Ra * psi + Lc * self.H(psi) - NI
        return f

    # Obtain h by using piecewise-linear interpolation
    def H(self, psi):
        b = psi / A
        i = 0
        for i in range(len(B) - 1):
            if b <= B[i + 1]:
                break
        x0 = B[i]
        x1 = B[i + 1]
        y0 = H[i]
        y1 = H[i + 1]
        m = (y1 - y0) / (x1 - x0)
        y = m * (b - x0) + y0
        return y

    # Calculate H prime
    def Hprime(self, psi):
        b = psi / A
        i = 0
        for i in range(len(B) - 1):
            if b <= B[i + 1]:
                break
        x0 = B[i]
        x1 = B[i + 1]
        y0 = H[i]
        y1 = H[i + 1]
        m = (y1 - y0) / (x1 - x0)
        return m

    # Function solves nonlinear equation using Successive Substitution method
    def successiveSubstitution(self, iguess, maxerror):
        count = 0
        psi = iguess
        fi = Ra * psi + Lc * self.H(iguess) - NI
        while (abs(self.f(psi) / fi) > maxerror):
            psi -= (Ra * psi + Lc * self.H(psi) - NI)*10**(-8)
            count += 1
        return psi, count
