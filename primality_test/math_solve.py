import random

def gcd(*nums):
    c = 0
    for x in nums:
        if c*x: 
            a,b = abs(x), c % abs(x)
            while b: a,b = b, a%b
            c = a
        else: c = abs(c+x)
    return c

def exgcd(a:int,b:int):
    if b == 0: return a,1,0
    c,x,y = exgcd(b,a%b)
    return c,y,x-(a//b)*y

def power_mod(a:int,b:int,N:int):
    assert b>=0, 'Power integer b should be NOT NEGATIVE'
    assert N>0, 'Modulo integer N should be POSITIVE'
    y = 1
    while True:
        if b % 2 > 0: y = y*a % N
        b //= 2
        if b == 0: break
        a = a*a % N
    return y

def Jacobi(a:int,n:int):
    assert n > 0 and n % 2 > 0, 'The second input integer should be POSITIVE and ODD'
    a %= n
    t = 1
    while a > 0:
        flip = n % 8 in [3,5]
        while a % 2 == 0:
            a //= 2
            if flip: t = -t
        if a % 4 == 3 and n % 4 == 3: t = -t
        a,n = n % a, a
    if n == 1: return t
    return 0

def Lucas_sequence(n:int,P:int,Q:int,N:int):
    assert n >= 0, 'Sequence number n should be NO NEGATIVE'
    assert N > 0, 'Modulo number N should be POSITIVE'
    if n == 0: return 0,2
    if n == 1: return 1,P
    ops = []
    while n > 1:
        ops.append(n%2)
        if n % 2 > 0: n -= 1
        else: n //= 2
    ops.reverse()
    D,U,V = P*P-4*Q,1,P
    for op in ops:
        if op: U,V = P*U+V,D*U+P*V
        else: U,V = 2*U*V,D*U*U+V*V
        U = (U+U%2*N)//2 % N
        V = (V+V%2*N)//2 % N
    return U,V

def powrt_int(n:int,k:int):
    if n == 1: return 1
    a,b = 1,n
    while b - a > 1:
        c = (a+b)//2
        y = c**k
        if y == n: return c
        if y > n: b = c
        else: a = c
    return 0

def trial_div_factorize(N:int):
    assert N > 0, 'Factorized integer should be POSITIVE.'
    if N == 1: return set()
    factors = set()
    y,m = N,2
    while y > 1 and m * m <= y:
        if y % m == 0:
            factors.add(m)
            while y % m == 0: y //= m
        m += 1 + m%2
    if y > 1: factors.add(y)
    return factors

def Pollard_rho(N:int,test_time=20,loop_time=2**16):
    assert N > 1, 'Factorized integer should be NO LESS THAN 2'
    for _ in range(test_time):
        c = random.randint(1,N-1)
        x1 = random.randint(0,N-1)
        x2 = (x1*x1+c) % N
        step = 0
        while x1 != x2 and step < loop_time:
            d = gcd(x1-x2,N)
            if d > 1:
                return d
            x1 = (x1*x1+c) % N
            x2 = (x2*x2+c) % N
            x2 = (x2*x2+c) % N
            step += 1
    return 1

def FR_factorize(N:int,test_time=10,balance=True):         # 默认尝试次数为10，变量balance表示是否希望F和R尽可能接近
    assert N > 1, 'Factorized integer should be NO LESS THAN 2'
    for _ in range(test_time):
        F = Pollard_rho(N)            # 先用Pollard's Rho算法初步找到一个因数
        if F == 1: break              # 如果找到的因数是1，说明可能无法分解或非常难分解，后面的尝试意义不大，直接结束
        R = N//F
        if F < R: F,R = R,F
        Q = gcd(F,R)
        while Q > 1:                # 要求F和R互素，如果不互素，就把公因数从F里拿出来乘到R上
            F //= Q
            R *= Q
            Q = gcd(F,R)
        if F < R: F,R = R,F        # 排个序，使F>R
        if R == 1: continue        # 如果得到的结果是 F=N, R=1，则本轮尝试无效，继续尝试
        if balance:                # 若balance为真，意味着希望F和R尽量接近，这样做是为了尽量减少后续分解F的计算量
            while R**3 < N:        # 当 R**3 < N 时，认为F与R相差较大
                F,Q = FR_factorize(F,test_time=test_time,balance=False)    # 继续分解F，这时不需要balance了，否则可能无限递归
                if Q == 1: break        # 如果F只能分解出F*1，则说明F可能已经无法分解了，不再尝试
                R *= Q                  # 若F能分解出非平凡因数，则将该因数拿出来乘到R上
        return (F,R) if F>R else (R,F)    # 由大到小返回F和R
    return N,1                            # 尝试次数到达上限，仍不能有效分解，则可能无法分解

def prob_prime_factorize(N:int,lim=2**50):
    if N <= lim: return trial_div_factorize(N),set()
    A = Pollard_rho(N)
    if A == 1: return set(),{N}
    B = N//A
    q1,q2 = prob_prime_factorize(A)
    q3,q4 = prob_prime_factorize(B)
    return q1|q3, q2|q4

def solve_quadratic_congruence(N,m):
    test_time = 100
    while True:
        test_time -= 1
        a = random.randint(0,N-1)
        d = (a*a-m) % N
        if Jacobi(d,N) == -1: break
        if test_time == 0: return {'info':'Jacobi(%d,N) == 1 but no solution for x^2 == %d mod N.'%(m,m)}
    n,u,v = N//2+1,a,1
    x,y = 1,0
    while True:
        if n % 2 > 0: x,y = (u*x+v*y*d) % N, (v*x+u*y) % N
        n //= 2
        if n == 0: break
        u,v = (u*u+v*v*d) % N, 2*u*v % N
    if y != 0 or (x*x-m) % N > 0: return {'info':'Jacobi(%d,N) == 1 but no solution for x^2 == %d mod N.'%(m,m)}
    return {x,N-x}



def Cornacchia_4N(N,D):
    xs = solve_quadratic_congruence(N,-D)
    if type(xs) == dict: return xs
    uvs = set()
    for x in xs:
        r,u = x,N%x
        while u*u >= N: r,u = u,r%u
        if (N-u*u) % D == 0:
            v2 = (N-u*u)//D
            v = powrt_int(v2,2)
            if v > 0: uvs.add((2*u,2*v))
            if D == 1 and v: uvs.add((2*v,2*u))
        if (-D) % 4 > 1: continue
        _,inv_4,inv_N = exgcd(4,N)
        ys = (x*4*inv_4+((-D)%4)*N*inv_N) % (4*N)
        ys = {ys,4*N-ys}
        for y in ys:
            r,u = y,(4*N)%y
            while u*u >= 4*N: r,u = u,r%u
            if (4*N-u*u) % D == 0:
                v2 = (4*N-u*u)//D
                v = powrt_int(v2,2)
                if v > 0: uvs.add((u,v))
                if D == 1 and v: uvs.add((v,u))
    return uvs

def V2_sequence(n,P1,P2,Q,N):
    if n == 0: return 2,0
    if n == 1: return 0,1
    if n == 2: return (-P2-2*Q)%N, P1%N
    inv_P2 = exgcd(P2,N)[1]
    # 定义递推关系式
    def double(x,y,a): return (x*x - P2*y*y - 2*power_mod(Q,a//2,N)) % N, (2*x + P1*y)*y % N
    def merge(x,y,r,s):
        v1 = (-Q*x-r)*inv_P2 % N
        return (Q*y+s-P1*v1) % N, v1
    def up(x,y,r,s):
        v0,v1 = merge(x,y,r,s)
        return (-P2*s-Q*v0) % N, (r+P1*s-Q*v1) % N
    # 计算递推路径 path_0,path_1
    path_0,path_1 = [],[n]
    while n > 2:
        if n % 2 > 0: n -= 1
        elif not path_0 or path_0[-1][0] == n-2: n //= 2
        if not path_0 and path_1[-1] % 2 == 0: path_1.append(n)
        elif path_0 and path_0[-1][1] - path_0[-1][0] == 2 : path_0.append((n-1,n))
        else: path_0.append((n-2,n))
    path_0.reverse()
    path_1.reverse()
    # 沿递推路径path_0进行递推
    for a,b in path_0:
        if b == 2: V0a,V1a,V0b,V1b = 2-2*a, a, (-P2-2*Q)%N, P1%N
        elif b - a == 2:
            V0a,V1a = double(V0a,V1a,a)
            V0b,V1b = double(V0b,V1b,b)
        elif a % 2 > 0: V0a,V1a = merge(V0a,V1a,V0b,V1b)
        else:
            x,y = up(V0a,V1a,V0b,V1b)
            V0a,V1a,V0b,V1b = V0b,V1b,x,y
    # 沿递推路径path_1进行递推
    for a in path_1:
        if a == 2: V0,V1 = (-P2-2*Q) % N, P1 % N
        elif a % 2 > 0: V0,V1 = up(V0a,V1a,V0b,V1b)
        else: V0,V1 = double(V0,V1,a)
    return V0,V1

def V3_sequence(n,P1,P2,P3,Q,N):
    if n == 0: return 2,0,0
    if n == 1: return 0,1,0
    if n == 2: return (-2*Q)%N, 0, 1
    inv_P3 = exgcd(P3,N)[1]
    # 定义递推关系式
    def double(x,y,z,a): return (x*x + 2*P3*y*z + P1*P3*z*z - 2*power_mod(Q,a>>1,N)) % N, (2*x*y - 2*P2*y*z + (P3-P1*P2)*z*z) % N, (y*y + 2*x*z + 2*P1*y*z + (P1*P1-P2)*z*z) % N
    def merge(x,y,z,r,s,t):
        v2 = (Q*x+r)*inv_P3 % N
        return (Q*y+s+P2*v2) % N, (Q*z+t-P1*v2) % N, v2
    def up(x,y,z,r,s,t):
        v0,v1,v2 = merge(x,y,z,r,s,t)
        return (P3*t-Q*v0) % N, (r-P2*t-Q*v1) % N, (s+P1*t-Q*v2) % N
    # 计算递推路径 path_0,path_1
    path_0,path_1 = [],[n]
    while n > 2:
        if n % 2 > 0: n -= 1
        elif not path_0 or path_0[-1][0] == n-2: n >>= 1
        if not path_0 and path_1[-1] % 2 == 0: path_1.append(n)
        elif path_0 and path_0[-1][1] - path_0[-1][0] == 2 : path_0.append((n-1,n))
        else: path_0.append((n-2,n))
    path_0.reverse()
    path_1.reverse()
    # 沿递推路径path_0进行递推
    for a,b in path_0:
        if b == 2: V0a,V1a,V2a,V0b,V1b,V2b = 2-2*a, a, 0, (-2*Q)%N, 0, 1
        elif b - a == 2:
            V0a,V1a,V2a = double(V0a,V1a,V2a,a)
            V0b,V1b,V2b = double(V0b,V1b,V2b,b)
        elif a % 2 > 0: V0a,V1a,V2a = merge(V0a,V1a,V2a,V0b,V1b,V2b)
        else:
            x,y,z = up(V0a,V1a,V2a,V0b,V1b,V2b)
            V0a,V1a,V2a,V0b,V1b,V2b = V0b,V1b,V2b,x,y,z
    # 沿递推路径path_0进行递推
    for a in path_1:
        if a == 2: V0,V1,V2 = (-2*Q)%N, 0, 1
        elif a % 2 > 0: V0,V1,V2 = up(V0a,V1a,V2a,V0b,V1b,V2b)
        else: V0,V1,V2 = double(V0,V1,V2,a)
    return V0,V1,V2

def primitive_root(q:int):
    ps = trial_div_factorize(q-1)
    for g in range(2,q):
        is_pt = True
        for p in ps:
            if power_mod(g,(q-1)//p,q) == 1:
                is_pt = False
                break
        if is_pt: return g

# 计算矩阵行列式，在求多项式判别式时会用到
# 这里采用最笨的递归法，对于小矩阵来说够用了
def det(M):
    d = len(M)
    if d == 1: return M[0][0]
    D = i = 0
    r = 1
    while i < d:
        if M[0][i]:
            M_ = [m[:i] + m[i+1:] for m in M[1:]]
            D += r * M[0][i] * det(M_)
        i += 1
        r = -r
    return D

def Euler_phi(n):
    S = R = n
    m = 2
    while m * m <= R:
        if R % m == 0: S = S*(m-1)//m
        while R % m == 0: R //= m
        #print(m,R,S)
        m += 1 + m%2
    if R > 1: S = S*(R-1)//R
    return S