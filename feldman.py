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
    return DSA.generate(1024)

class Diler:
    def __init__(self):
        self._p = None
        self._q = None
        self._n = 5
        self._t = 3
        self._g = None
        self._k = randrange(100, self._p - 100)
        self._coeff = self.get_coeff()
        print('current secret: ', self._k)

    def set_params(self):
        keys = get_params()
        self._p = keys.p#90788617789067741660093735441816283768126454165447820914348986349362240012246272190503132710290350012724240057975518271126302221327143516242885041310935673700665712867659263485029894251155508346242612663272171256883981343369810628426686764826164140935803164377306581161297913701577432630419189938806935481267
        self._q = keys.q#11858713048509395322704712338910337791253467170330157609738243194501589467866347742859026652173247929744372823851475690076568633010816669049447174858448081
        self._g = keys.g#13495691096980511632596807502578684692636842856060619187984980202240173176520436277561555543447166164836371310304905899110868343032793007222209286409424594712307915797199849718973347960967229169667233906242901259266594119591586764827309527398369094753391559746710037551737980047139900412522248101415392617175

    def get_coeff(self):
        self.set_params()
        coeff = [randrange(self._p) for _ in range(self._t - 1)]
        coeff.append(self._k)
        #print(coeff)
        return coeff

    def get_polynom(self, x):
        beta = 0
        for coeff_ind, coeff_val in enumerate(self._coeff[::-1]):
            beta += int(gmpy2.powmod(x, coeff_ind, self._p)) * coeff_val

        return beta % self._q

    def set_prooves(self):
        self._prooves = [int(gmpy2.powmod(self._g, coeff_val, self._p)) % self._p for _, coeff_val in enumerate(self._coeff[::-1])]

    def get_prooves(self):
        self.set_prooves()
        return self._prooves

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
                    A[i][j] = int(gmpy2.powmod(parts[i][0], (j+1), self._p))
        #A = GF([[parts[0][0], parts[0][0]**2, 1], [parts[1][0], parts[1][0]**2, 1], [parts[2][0], parts[2][0]**2, 1]])
        b_list = [parts[i][1] for i in range(self._t)]
        b = GF(b_list)
        print('A: ', A)
        print('b: ', b)
        secret = np.linalg.solve(A, b)
        print('current secret:', secret[-1])

    def check_prooves(self, part, prooves):
        #check = [gmpy2.powmod(self._g, gmpy2.powmod(parts[0][0], i, self._p),self._p) for coeff_ind, coeff_val in enumerate(self._coeff[::-1]) )]
        z = int(gmpy2.powmod(self._g, part[1], self._p))
        #check = mult([int(gmpy2.powmod(proove, gmpy2.powmod(part[0], i, self._p), self._p)) for i, proove in
        #         enumerate(prooves)])
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

if __name__ == "__main__":
    D = Diler()
    print('p: ', D._p)
    print('q: ', D._q)
    print('g: ', D._g)
    parts = []
    parts = D.generate_shares()
    #print(parts)
    prooves = D.get_prooves()
    D.check_prooves(parts[2], prooves)
    check_parts = [parts[0], parts[2], parts[4]]
    print(D.lagr_recon_secret(check_parts))
    print('p size: ', len(str(D._q)))
    #D.scnd_recon_secret(check_parts)
    '''k = get_params()
    print(k.p)
    print(k.q)
    print(k.g)'''
    print((D._p - 1) % D._q)
    print(gmpy2.powmod(D._g,D._q,D._p))