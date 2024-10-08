import math
import numpy
import random
from .math_solve import gcd,exgcd,Jacobi,solve_quadratic_congruence
from .polynomial import Polynomial,poly_pow_mod,gcd_poly

def prime_list_for_schoof(N:int) -> list:
    if N == 5: return [7]
    n = int(math.log(4*math.sqrt(N)))+100
    primality = [True]*(n+1)
    primes = []
    T = m = 2
    while m * m <= N:
        for k in range(2,n//m+1): primality[k*m] = False
        m += 1
        while not primality[m]: m +=1
        primes.append(m)
        T *= m
        if T*T > 16*N: return primes
    m += 1
    while m <= n:
        if primality[m]:
            primes.append(m)
            T *= m
            if T*T > 16*N: return primes
        m += 1

def class_group(D):
    matrices = set()
    B = 0
    while 3*B*B < D:
        A = 2*B
        while A*A <= D+B*B:
            if A > 0 and (D+B*B) % A == 0:
                C = (D+B*B) // A
                if gcd(A,2*B,C) == 1:
                    matrices.add((A,B,C))
                    if 0 < 2*B < A < C: matrices.add((A,-B,C))
            A += 1
        B += 1
    return matrices

def reduced_class_polynomial(D):
    def F(z): return 1 - z - z*z + z**5 + z**7 - z**12 - z**15 #1-z*(1+z*(1-z**3*(1+z*z*(1-z**5*(1+z**3)))))
    classgroup = class_group(D)
    degree = len(classgroup)
    HD = numpy.zeros(degree+1)*0j
    HD[0] = 1
    G = 1 if D%3>0 else 3
    I = 3 if D%8 in [1,2,6,7] else (6 if D%8==5 else (0 if D%3>0 else 2))
    K = 4 if D%8==5 else (1 if D%4==3 else 2)
    for A,B,C in classgroup:
        L = (A-C+A*C*C*(5 if D%8==3 else -1)) if A%2==0 else ((A-C+A*A*C) if C%2>0 or D%8==5 else (A+2*C-A*C*C))
        M = (-1) if A%2>0 and A*A//8%2>0 or A%2==0 and C*C//8%2>0 else 1
        N = (-M) if D%8==3 and A*C%2==0 else (M if D%8 in [1,2,6] or D%8==7 and A*C%2>0 else 1)
        theta = numpy.exp((-numpy.sqrt(D)+B*1j)*numpy.pi/A)
        if A*C % 2 > 0: f = theta**(-1/24) * F(-theta) / F(theta**2)
        elif A % 2 > 0: f = theta**(-1/24) * F(theta) / F(theta**2)
        else: f = numpy.sqrt(2) * theta**(1/12) * F(theta**4) / F(theta**2)
        CI = (N*numpy.exp(1j*numpy.pi*K*B*L/24)*2**(-I/6)*f**K)**G
        HD[1:] = HD[:-1]
        HD[0] = 0
        HD[:-1] -= CI*HD[1:]
    coef = {}
    for i in range(degree+1):
        c = int(HD[i].real)
        if abs(c-HD[i]) < 1e-5: coef[i] = c
        elif abs(c+1-HD[i]) < 1e-5: coef[i] = c+1
        elif abs(c-1-HD[i]) < 1e-5: coef[i] = c-1
        else: return {'info':'Non-Integer coefficient occurs in reduced class polynomial ...'}
    return Polynomial(coef)

def complex_multiplication(N,D,RCP,u_odd):
    if D == 1: return 1,0
    if D == 2: return (-30)%N, 56%N
    if D == 3: return 0,1
    if D == 7: return (-35)%N, 98%N
    if D == 11: return (-264)%N, 1694%N
    if D == 19: return (-152)%N, 722%N
    if D == 43: return (-3440)%N, 77658%N
    if D == 67: return (-29480)%N, 1948226%N
    if D == 163: return (-8697680)%N, 9873093538%N
    if u_odd:
        g = RCP.find_factor(3)
        Z = (Polynomial({24:-1},N) if D%3>0 else Polynomial({8:-256},N)) % g
        a1 = (-3) * (Z+64) * (Z+256) % g
        b1 = 2 * (Z+64)**2 * (Z-512) % g
        a3 = a1**3 % g
        b2 = b1*b1 % g
        for k in range(max(a3.degree,b2.degree)+1):
            if k in a3.coef.keys() and k in b2.coef.keys():
                sigma,tau = a3.coef[k],b2.coef[k]
                break
        a = sigma*tau % N
        b = a*tau % N
    else:
        s = -(RCP.find_factor(1).coef[0])
        G = gcd(D,3)
        I = 3 if D%8 in [1,2,6,7] else (6 if D%8==5 else (0 if D%3>0 else 2))
        K = 4 if D%8==5 else (1 if D%4==3 else 2)
        Z = (-1)**(D&1) * 2**(4*I//K) * s**(24//(G*K)) % N
        a = (-3) * (Z+64) * (Z+16) % N
        b = 2 * (Z+64)**2 * (Z-8) % N
    return a,b

class Elliptic_curve():
    def __init__(self,N,a,b):
        self.N = N
        self.a = a
        self.b = b
    def point_on_curve(self,x,y):
        return (y*y - x**3 - self.a*x - self.b) % self.N == 0
    def random_point(self,zero_point=True):          # if zero_point is true, points on x-axis are acceptable
        while True:
            x = random.randint(0,self.N-1)
            y2 = (x**3 + self.a*x + self.b) % self.N
            if y2 == 0:
                if zero_point: return (x,0)
                continue
            jacobi = Jacobi(y2,self.N)
            if jacobi == -1: continue
            if jacobi == 0: return {'info':'Factor found: %d'%gcd(y2,self.N)}        # nontrivial factor of N is found
            ys = solve_quadratic_congruence(self.N,y2)
            if type(ys) == dict: return ys                                # information that N is composite
            return (x,list(ys)[random.randint(0,1)])
    def plus(self,P,Q):
        if P == (0,): return Q    # tuple (0,) represents infinity point
        if Q == (0,): return P
        if P[0] != Q[0]:
            c,inv_dx,_ = exgcd(P[0]-Q[0],self.N)
            if c > 1: return {'info':'Factor found: %d'%c}   # nontrivial factor of N is found
            k = (P[1]-Q[1]) * inv_dx % self.N
            x = (k*k-P[0]-Q[0]) % self.N
            y = (-k*(x-P[0])-P[1]) % self.N
            if not self.point_on_curve(x,y): return {'info':'Group operation failed: point not on the curve'}
            return (x,y)
        if (P[1] + Q[1]) % self.N == 0: return (0,)
        if P[1] == Q[1]:
            c,inv_2y,_ = exgcd(2*P[1],self.N)
            if c > 1: return {'info':'Factor found: %d'%c}   # nontrivial factor of N is found
            k = ((3*P[0]*P[0]+self.a)*inv_2y) % self.N
            x = (k*k-2*P[0]) % self.N
            y = (-k*(x-P[0])-P[1]) % self.N
            if not self.point_on_curve(x,y): return {'info':'Group operation failed: point not on the curve'}
            return (x,y)
        return {'info':'Not a group: x1==x2 but y1!=y2 and y1+y2!=0'}
    def times(self,m,P):     # an integer m multiples a point P
        Q = (0,)
        while True:
            if m % 2 > 0: Q = self.plus(Q,P)
            if type(Q) == dict: return Q    # information that N is composite
            m //= 2
            if m == 0: break
            P = self.plus(P,P)
            if type(P) == dict: return Q    # same as above
        return Q
    def schoof(self):                 # Schoof's algorithm
        def const(c): return Polynomial(c,self.N)
        primes = prime_list_for_schoof(self.N)
        x = Polynomial(modulo=self.N)
        E = x**3 + self.a*x + self.b
        xN = poly_pow_mod(x,self.N,E)
        if type(xN) == dict: return xN        # information that N is composite
        pcd = gcd_poly(E,xN-x)
        if type(pcd) == dict: return pcd
        t = int(not pcd.degree)
        M = 2
        f = [const(0), const(1), const(2), 3*x**4 + (6*self.a%self.N)*x*x + (12*self.b%self.N)*x - self.a**2%self.N,
             4*x**6 + (20*self.a%self.N)*x**4 + (80*self.b%self.N)*x**3 - (20*self.a**2%self.N)*x*x - (16*self.a*self.b%self.N)*x - (4*self.a**3+32*self.b**2)%self.N]
        yj0,yj1 = E*E,1
        inv_2 = exgcd(2,self.N)[1]
        for j in range(2,primes[-1]//2+2):
            f.append(f[j]**3*f[j+2]*yj0-f[j-1]*f[j+1]**3*yj1)
            f.append(f[j+1]*(f[j]**2*f[j+3]-f[j-1]*f[j+2]**2)*inv_2)
            yj0,yj1 = yj1,yj0
        f.append(const(-1))
        for m in primes:
            xN = poly_pow_mod(x,self.N,f[m])
            xN2 = poly_pow_mod(x,self.N**2,f[m])
            yN = poly_pow_mod(E,self.N//2,f[m])
            if type(xN) == dict: return xN
            if type(xN2) == dict: return xN2
            if type(yN) == dict: return yN
            s = self.N % m
            ys0,ys1 = (1,E) if s%2>0 else (E,1)
            g1 = (xN2-x)*f[s]**2*ys0 + f[s-1]*f[s+1]*ys1
            pcd = gcd_poly(g1,f[m])
            if type(pcd) == dict: return pcd
            if pcd.degree:
                if Jacobi(s,m) == -1: tm = 0
                else:
                    for w in range(1,m//2+1):
                        if (w*w-s) % m == 0: break
                    yw0,yw1 = (1,E) if w%2>0 else (E,1)
                    g2 = (xN-x)*f[w]**2*yw0 + f[w-1]*f[w+1]*yw1
                    pcd = gcd_poly(g2,f[m])
                    if type(pcd) == dict: return pcd
                    if not pcd.degree: tm = 0
                    else:
                        g3 = f[w-2]*f[w+1]**2 - f[w-1]**2*f[w+2] + 4*f[w]**3*yN*yw0**2
                        pcd = gcd_poly(g3,f[m])
                        if type(pcd) == dict: return pcd
                        tm = 2*w if pcd.degree>0 else -2*w
            else:
                yN2 = poly_pow_mod(E,self.N**2//2,f[m])
                if type(yN2) == dict: return yN2
                alpha = f[s-1]**2*f[s+2] - f[s-2]*f[s+1]**2 - 4*f[s]**3*yN2*ys0**2
                beta = -4*f[s]*g1
                h1 = (alpha**2*ys1-(xN2+xN+x)*beta**2*ys0)*f[s]**2 + beta**2*f[s-1]*f[s+1]*ys1
                h2 = beta**2 * f[s]**2 * ys0
                tau = 1
                while tau <= (m+1)//2:
                    if tau % 2 > 0: ft1,ft2 = poly_pow_mod(f[tau],2*self.N,f[m]),poly_pow_mod(f[tau-1]*f[tau+1]*E,self.N,f[m])
                    else: ft1,ft2 = poly_pow_mod(f[tau]**2*E,self.N,f[m]),poly_pow_mod(f[tau-1]*f[tau+1],self.N,f[m])
                    if type(ft1) == dict: return ft1
                    if type(ft2) == dict: return ft2
                    g4x = h1 * ft1 + h2 * ft2
                    pcd = gcd_poly(g4x,f[m])
                    if type(pcd) == dict: return pcd
                    if pcd.degree:
                        h3 = ((2*xN2+x)*alpha*beta**2*ys0-alpha**3*ys1-beta**3*yN2*ys0**2)*f[s]**2 - alpha*beta**2*f[s-1]*f[s+1]*ys1
                        h4 = beta**3 * f[s]**2 * ys0**2
                        ft3,ft4 = poly_pow_mod(f[tau],3*self.N,f[m]),poly_pow_mod(f[tau-1]**2*f[tau+2]-f[tau-2]*f[tau+1]**2,self.N,f[m])
                        if type(ft3) == dict: return ft3
                        if type(ft4) == dict: return ft4
                        if tau % 2 > 0: g4y = 4*h3*ft3 - h4*ft4*yN
                        else: g4y = 4*h3*ft3*yN**3*E*E - h4*ft4
                        pcd = gcd_poly(g4y,f[m])
                        if type(pcd) == dict: return pcd
                        tm = tau if pcd.degree>0 else -tau
                        break
                    tau += 1
                if tau > (m+1)//2: return {'info':'Frobenius endomorphism failed.'}   # information that N is composite
            _,inv_M,inv_m = exgcd(M,m)
            t = (tm*M*inv_M+t*m*inv_m) % (m*M)
            M *= m
        if t >= M//2: t -= M
        if t*t > 4*self.N: return {'info':'Hasse theorem failed.'}      # information that N is composite
        return self.N+1-t
