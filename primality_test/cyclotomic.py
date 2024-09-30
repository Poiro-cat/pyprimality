from .math_solve import power_mod,exgcd

class Cyclotomic():
    def __init__(self,p:int,k=1,coef={}):
        self.p = p
        self.k = k
        self.num = p**k
        self.coef = [0] * self.num     # 储存各个次项的整系数
        if type(coef) == dict:
            for i in coef: self.coef[i%self.num] = coef[i]
        if type(coef) == list:
            assert self.num == len(coef), 'Length of coef should be equal to p^k.'
            self.coef = coef
    def copy(self):
        return Cyclotomic(self.p,self.k,self.coef.copy())
    def __add__(self,a):        # 加法
        z = self.copy()
        if type(a) == int: z.coef[0] += a
        if type(a) == Cyclotomic:
            assert z.num == a.num, 'Dimension of two instances should be the same.'
            for i in range(self.num): z.coef[i] += a.coef[i]
        return z
    def __radd__(self,a):        # 加法
        return self + a
    def __neg__(self):           # 负号
        z = self.copy()
        for i in range(self.num): z.coef[i] *= -1
        return z
    def __sub__(self,a):         # 减法
        return self + (-a)
    def __rsub__(self,a):        # 减法
        return (-self) + a
    def __mul__(self,a):         # 乘法
        z = Cyclotomic(self.p,self.k)
        if type(a) == int:
            for i in range(self.num): z.coef[i] = self.coef[i] * a
        if type(a) == Cyclotomic:
            assert z.num == a.num, 'Dimension of two instances should be the same.'
            for i in range(self.num):
                for j in range(a.num):
                    r = (i+j) % self.num
                    z.coef[r] += self.coef[i]*a.coef[j]
        return z
    def __rmul__(self,a):         # 乘法
        return self * a
    def __mod__(self,modulo):       # 模除法，令每个整系数取模
        z = self.copy()
        for i in range(self.num): z.coef[i] %= modulo
        return z
    def __pow__(self,n):         # 幂运算
        z = 1
        y = self.copy()
        while True:
            if n % 2 > 0: z *= y
            n //= 2
            if n == 0: return z
            y *= y
    def sigma(self,x:int):        # 算符sigma_x的运算
        z = Cyclotomic(self.p,self.k)
        for i in range(self.num):
            j = (i*x) % self.num
            z.coef[j] += self.coef[i]
        return z
    def ZG(self,G:dict,modulo=0):     # 算符sigma_x线性组合的指数运算
        z = 1
        for i in G:
            y = self.sigma(i)
            z *= power_mod(y,G[i],modulo) if modulo>0 else y**G[i]
            if modulo: z %= modulo
        return z
    def simplify(self):           # 利用复数单位圆的性质化简
        z = self.copy()
        t = z.num // z.p
        for i in range(z.num-t,z.num):
            for j in range(i+t-z.num,i,t): z.coef[j] -= z.coef[i]
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
            if k > 1: st = st + '^' + str(k)
        st = st[3:] if st[1]=='+' else st[1:]
        return st

# 计算 j_pq，当 p^k > 2 时
def calculate_jpq(p,k,q,f):
    j = Cyclotomic(p,k)
    for x in range(1,q-1): j.coef[(x+f[x])%j.num] += 1
    return j

# 计算 j_2q(*) 和 j_2q(#)，当 p == 2 且 k > 2 时
def calculate_j2q(k,q,j,f):
    j1,j2 = Cyclotomic(2,k),Cyclotomic(2,k)
    for x in range(1,q-1):
        j1.coef[(2*x+f[x])%j1.num] += 1
        j2.coef[(3*x+f[x]<<k-3)%j2.num] += 1
    return j1 * j, j2 * j2

# 计算 J_pq，当 p > 2 时
def calculate_Jpq(p,k,vk,j,N):
    M = [i//(p-1)*p + i%(p-1) + 1 for i in range((p-1)*p**(k-1))]
    theta,alpha = {},{}
    for x in M:
        y = exgcd(x,p**k)[1] % p**k
        theta[y] = x
        alpha[y] = vk*x // p**k
    return j.ZG(theta,N), j.ZG(alpha,N)

# 计算 J_2q，当 p == 2 且 k > 2 时
def calculate_J2q(k,vk,j_ast,j_tag,N):
    theta,alpha = {},{}
    for x in range(1,2**k,2):
        if x % 8 > 4: continue
        y = exgcd(x,1<<k)[1] % 2**k
        theta[y] = x
        alpha[y] = vk*x // 2**k
    return j_ast.ZG(theta,N), j_ast.ZG(alpha,N)*(j_tag if vk%8>4 else 1)

# 在步骤4.3和5.4中寻找满足条件的整数h
def find_h(p,k,z,N):
    h,y = 0,Cyclotomic(p,k,{0:1})
    zeta = Cyclotomic(p,k,{1:1})
    while h < p**k:
        z1 = (z-y).simplify() % N
        if not any(z1.coef): break
        h += 1
        y *= zeta
    return h