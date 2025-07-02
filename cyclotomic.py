from .utils import power_mod, exgcd

class Cyclotomic():
    def __init__(self, p:int, k:int=1, coef={}):
        self.p = p
        self.k = k
        self.num = p**k
        self.coef = [0] * self.num
        if type(coef) == dict:
            for i in coef: self.coef[i%self.num] = coef[i]
        elif type(coef) == list:
            assert self.num == len(coef),  'Length of coef should be equal to p^k.'
            self.coef = coef
        else: raise TypeError('Input coef must be a list or a dictionary.')
    def copy(self):
        return Cyclotomic(self.p, self.k, self.coef.copy())
    def __add__(self, a):
        z = self.copy()
        if type(a) == int: z.coef[0] += a
        elif type(a) == Cyclotomic:
            assert z.num == a.num,  'Dimension of two instances should be the same.'
            for i in range(self.num): z.coef[i] += a.coef[i]
        else: raise TypeError('Operated instances must be Integer or Cyclotomic.')
        return z
    def __radd__(self, a:int):
        return self + a
    def __neg__(self):
        z = self.copy()
        for i in range(self.num): z.coef[i] *= -1
        return z
    def __sub__(self, a):
        return self + (-a)
    def __rsub__(self, a):
        return (-self) + a
    def __mul__(self, a):
        z = Cyclotomic(self.p, self.k)
        if type(a) == int:
            for i in range(self.num): z.coef[i] = self.coef[i] * a
        elif type(a) == Cyclotomic:
            assert z.num == a.num,  'Dimension of two instances should be the same.'
            for i in range(self.num):
                for j in range(a.num):
                    r = (i+j) % self.num
                    z.coef[r] += self.coef[i]*a.coef[j]
        else: raise TypeError('Operated instances must be Integer or Cyclotomic.')
        return z
    def __rmul__(self, a:int):
        return self * a
    def __mod__(self, modulo:int):
        z = self.copy()
        for i in range(self.num): z.coef[i] %= modulo
        return z
    def __pow__(self, n:int):
        z = 1
        y = self.copy()
        while True:
            if n % 2 > 0: z *= y
            n //= 2
            if n == 0: return z
            y *= y
    def sigma(self, x:int):        # computation of operator sigma_x
        z = Cyclotomic(self.p, self.k)
        for i in range(self.num):
            j = (i*x) % self.num
            z.coef[j] += self.coef[i]
        return z
    def ZG(self, G:dict, modulo:int=0):     # computation of linear combination of sigma_x's
        z = 1
        for i in G:
            y = self.sigma(i)
            z *= power_mod(y, G[i], modulo) if modulo>0 else y**G[i]
            if modulo: z %= modulo
        return z
    def simplify(self):           # using the property of unit circle in the complex filed
        z = self.copy()
        t = z.num // z.p
        for i in range(z.num-t, z.num):
            for j in range(i+t-z.num, i, t): z.coef[j] -= z.coef[i]
            z.coef[i] = 0
        return z
    def __str__(self):
        if not any(self.coef): return '0'
        st = ''
        for k in range(self.num):
            c = self.coef[k]
            if not c: continue
            if c > 0: st = st + ' + '
            if c < 0: st = st + ' - '
            if not k or c*c > 1: st = st + str(c if c > 0 else -c)
            if c*c > 1 and k: st = st + '*'
            if k: st = st + 'Z'
            if k > 1: st = st + f'^{k}'
        st = st[3:] if st[1]=='+' else st[1:]
        return st

# compute j_pq, when p^k > 2
def calculate_jpq(p:int, k:int, q:int, f:list[int]):
    j = Cyclotomic(p, k)
    for x in range(1, q-1): j.coef[(x+f[x]) % j.num] += 1
    return j

# compute j_2q(*) and j_2q(#), when p == 2 and k > 2
def calculate_j2q(k:int, q:int, j:Cyclotomic, f:list[int]):
    j1, j2 = Cyclotomic(2, k), Cyclotomic(2, k)
    for x in range(1, q-1):
        j1.coef[(2*x+f[x]) % j1.num] += 1
        j2.coef[(3*x+f[x]) * 2**(k-3) % j2.num] += 1
    return j1 * j,  j2 * j2

# compute J_pq, when p > 2
def calculate_Jpq(p:int, k:int, vk:int, j:Cyclotomic, N:int):
    M = [i//(p-1)*p + i%(p-1) + 1 for i in range((p-1)*p**(k-1))]
    theta, alpha = {}, {}
    for x in M:
        y = exgcd(x, p**k)[1] % p**k
        theta[y] = x
        alpha[y] = vk*x // p**k
    return j.ZG(theta, N),  j.ZG(alpha, N)

# compute J_2q, when p == 2 and k > 2
def calculate_J2q(k:int, vk:int, j_ast:Cyclotomic, j_tag:Cyclotomic, N:int):
    theta, alpha = {}, {}
    for x in range(1, 2**k, 2):
        if x % 8 > 4: continue
        y = exgcd(x, 1<<k)[1] % 2**k
        theta[y] = x
        alpha[y] = vk*x // 2**k
    return j_ast.ZG(theta, N),  j_ast.ZG(alpha, N)*(j_tag if vk%8>4 else 1)

# find integer h in step 4.3 and 5.4 of APR-CL
def find_apr_h(p:int, k:int, z:Cyclotomic, N:int):
    h, y = 0, Cyclotomic(p, k, {0:1})
    zeta = Cyclotomic(p, k, {1:1})
    while h < p**k:
        z1 = (z-y).simplify() % N
        if not any(z1.coef): break
        h += 1
        y *= zeta
    return h
