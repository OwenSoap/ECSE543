# Author@Yi Zhu
# ID@260716006

class Integral:

    # Function calculate integral by using Gauss-Legendre method with uniform segments
    def gaussLegendreUni(self, f, n, a, b):
        n, a, b = float(n), float(a), float(b)
        sum = 0
        w = (b - a) / n
        h = [w] * int(n)
        for w in h:
            low = a
            a += w
            up = a
            sum += (up - low) * f((low + up) / 2)
        return sum

    # Function calculate integral by using Gauss-Legendre method with non-uniform segments
    def gaussLegendreNonUni(self, f, segs):
        sum = 0
        h = len(segs)
        for i in range(1, h):
            b = segs[i]
            a = segs[i - 1]
            sum += (b - a) * f((a + b) / 2)
        return sum
