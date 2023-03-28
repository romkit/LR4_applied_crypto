from Crypto.Random.random import randrange
from Crypto.Util.number import inverse, getStrongPrime, getPrime
from Crypto.PublicKey import DSA
from functools import reduce
import numpy as np
import galois
import gmpy2

def mult(lst):
    return reduce(lambda x, y: x * y, lst, 1)

def get_params():
    keys = DSA.generate(1024)
    d = randrange(100, keys.q)
    h = gmpy2.powmod(keys.g, d, keys.p)
    return keys, d, h


class Diler:
    def __init__(self):
        self._p = None
        self._q = None
        self._n = 5
        self._t = 3
        self._g = None
        self._k = None
        self._coeff_beta = list()
        self._coeff_gamma = list()
        self._d = None
        self._h = None
        self._gam = None #randrange(100, self._p - 100)
        #print('current secret: ', self._k)

    def set_params(self):
        keys, self._d, self._h = get_params()
        self._p = keys.p
        self._q = keys.q
        self._g = keys.g
        self._k = randrange(100, self._p - 100)
        self._coeff_beta = self.set_coeff(self._k)
        self._coeff_gamma = self.set_coeff(randrange(100))
        print('current secret: ', self._k)

    def set_coeff(self, first):
        coeff = [randrange(self._p) for _ in range(self._t - 1)]
        coeff.append(first)
        #print(coeff)
        return coeff

    def get_polynoms(self, x):
        #self.set_params()
        beta = 0
        for coeff_ind, coeff_val in enumerate(self._coeff_beta[::-1]):
            beta += int(gmpy2.powmod(x, coeff_ind, self._p)) * coeff_val

        gamma = 0
        for coeff_ind, coeff_val in enumerate(self._coeff_gamma[::-1]):
            gamma += int(gmpy2.powmod(x, coeff_ind, self._p)) * coeff_val


        return beta % self._q, gamma % self._q

    def get_prooves(self):
        eps = list()
        for i in range(1, self._t + 1):
            beta, gamma = self.get_polynoms(i)
            eps.append(int(gmpy2.powmod(self._g, beta, self._p) *
                  gmpy2.powmod(self._h, gamma, self._p)) % self._p)
        return eps

    def get_shares(self):
        parts = list()
        for i in range(1, self._n + 1):
            x = i
            parts.append((x, self.get_polynoms(x)))

        return parts

    def check_prooves(self, part, prooves):
        z = int(gmpy2.powmod(self._g, part[1][0], self._p) * gmpy2.powmod(self._h, part[1][1], self._p))
        check = list()
        for i, proove in enumerate(prooves):
            pre_tmp = gmpy2.powmod(part[0], i, self._p)
            tmp = int(gmpy2.powmod(proove, pre_tmp, self._p))
            check.append(tmp)
        check = mult(check)
        if z % self._p == check % self._p:
            print('true checker: ')
        else:
            print('false checker')

    def lagr_recon_secret(self, parts):
        summary = 0
        for j, part_j in enumerate(parts):
            x_j, y_j = part_j
            prod = 1
            for i, part_i in enumerate(parts):
                x_i, _ = part_i
                if i != j:
                    prod *= x_i * inverse(x_i - x_j, self._p)

            prod *= y_j[0]
            #print(prod)
            summary += prod

        if summary % self._p != self._k:
            print('false reset key')
        return summary % self._p



if __name__ == "__main__":
    D = Diler()
    D.set_params()
    print('p: ', D._p)
    print('q: ', D._q)
    print('g: ', D._g)
    parts = []
    parts = D.get_shares()
    print('parts: ', parts)
    prooves = D.get_prooves()
    #print('prooves: ', prooves)
    #check_parts = [parts[0], parts[2], parts[3]]
    #print(D.lagr_recon_secret(check_parts))
    D.check_prooves(parts[2], prooves)

