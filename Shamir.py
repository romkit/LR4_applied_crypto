# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
from Crypto.Random.random import randrange
from Crypto.Util.number import inverse, getStrongPrime, getPrime

import numpy as np
import galois
import gmpy2

class Diler:
    def __init__(self):
        self._p = getPrime(512)
        self._n = 5
        self._t = 3
        self._k = randrange(100, self._p - 100)
        self._coeff = self.get_coeff()
        print('current secret: ', self._k)

    def get_coeff(self):
        coeff = [randrange(self._p) for _ in range(self._t - 1)]
        coeff.append(self._k)
        #print(coeff)
        return coeff

    def get_polynom(self, x):
        point = 0
        for coeff_ind, coeff_val in enumerate(self._coeff[::-1]):
            point += int(gmpy2.powmod(x, coeff_ind, self._p)) * coeff_val
        return point % self._p

    def generate_shares(self):
        parts = list()

        for i in range(1, self._n + 1):
            x = i
            parts.append((x, self.get_polynom(x)))

        return parts

    def lagr_recon_secret(self, parts):
        summary = 0
        for j, part_j in enumerate(parts):
            x_j, y_j = part_j
            prod = 1
            for i, part_i in enumerate(parts):
                x_i, _ = part_i
                if i != j:
                    prod *= x_i * inverse(x_i - x_j, self._p)

            prod *= y_j
            #print(prod)
            summary += prod
        return summary % self._p

    def scnd_recon_secret(self, parts):
        GF = galois.GF(self._p)
        A = GF.Identity(self._t)

        for i in range(self._t):
            for j in range(self._t):
                if j == self._t - 1:
                    A[i][j] = 1
                else:
                    A[i][j] = parts[i][0] ** (j+1)
        #A = GF([[parts[0][0], parts[0][0]**2, 1], [parts[1][0], parts[1][0]**2, 1], [parts[2][0], parts[2][0]**2, 1]])
        b_list = [parts[i][1] for i in range(self._t)]
        b = GF(b_list)
        #print(A)
        #print(b)
        secret = np.linalg.solve(A, b)
        print('current secret:', secret[-1])


if __name__ == "__main__":
    D = Diler()
    print(D._p)
    parts = []
    parts = D.generate_shares()
    print(parts)
    check_parts = [parts[0], parts[2], parts[3]]
    print(D.lagr_recon_secret(check_parts))
    #D.scnd_recon_secret(check_parts)