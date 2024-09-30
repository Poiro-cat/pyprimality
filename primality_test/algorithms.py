import math
import random

from .math_solve import (gcd,exgcd,power_mod,Jacobi,Lucas_sequence,powrt_int,
                         trial_div_factorize,Pollard_rho,FR_factorize,prob_prime_factorize,
                         Cornacchia_4N,V2_sequence,V3_sequence,Euler_phi)
from .elliptic_curve import reduced_class_polynomial,complex_multiplication,Elliptic_curve
from .cyclotomic import calculate_jpq,calculate_j2q,calculate_Jpq,calculate_J2q,find_h
from .apr_functions import trial_division_30,v,e,generate_f
from .polynomial import Polynomial,poly_pow_mod,gcd_poly,poly_nest


##### 0. Trial division #####

def Trial_division(N:int,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    m = 2
    while m * m <= N:
        if N % m == 0:
            if pr: print('Factor found: %d'%m)
            if pr: print('COMPOSITE')
            return False
        m += 1
    if pr: print('No factor found.')
    if pr: print('PRIME')
    return True

def Trial_division_30(N:int,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    factor = 1
    for m in [2,3,5]:
        if m * m > N and factor == 1: break
        if N % m == 0: factor = m
    if factor == 1 and m * m < N:
        A,i,m = 0,1,7
        B = [1,7,11,13,17,19,23,29]
        while m * m <= N:
            if N % m == 0:
                factor = m
                break
            i = (i+1) % 8
            if i == 0: A += 30
            m = A + B[i]
    if factor > 1:
        if pr: print('Factor found: %d'%factor)
        if pr: print('COMPOSITE')
        return False
    if pr: print('No factor found.')
    if pr: print('PRIME')
    return True

def Trial_division_lim(N:int,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    primes = []
    m,A,B,k,i,j = 2,2,[1],1,1,0
    while True:
        #print('Try: m =',m)
        if m * m > N:
            if pr: print('No factor found.')
            if pr: print('PRIME')
            return True
        primality = True
        for p in primes:
            if p * p > m: break
            if m % p == 0: primality = False
        if primality:
            primes.append(m)
            if N % m == 0:
                if pr: print('Factor found: %d'%m)
                if pr: print('COMPOSITE')
                return False
        m = A*k + B[j]
        j += 1
        if j >= len(B):
            k += 1
            j = 0
        if m > 3 and k >= primes[i]:
            c = primes[i]
            B_ = []
            for k in range(c):
                for b in B:
                    if (A*k+b) % c > 0: B_.append(A*k+b)
            A *= c
            i += 1
            k = 1
            B,j = B_,0


##### 1. Fermat #####

def Fermat(N:int,test_time=100,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    for _ in range(test_time):
        a = random.randint(2,N-2)
        if gcd(a,N) > 1:
            if pr: print('Factor found: %d'%gcd(a,N))
            if pr: print('COMPOSITE')
            return False
        if power_mod(a,N-1,N) != 1:
            if pr: print('Fermat test failed with a = %d'%a)
            if pr: print('COMPOSITE')
            return False
    if pr: print('Fermat test passed with %d repeated tests.'%test_time)
    if pr: print('PRIME (probable)')
    return True


##### 2. Miller-Rabin #####

def Miller_Rabin(N:int,test_time=100,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    for _ in range(test_time):
        a = random.randint(2,N-2)
        if gcd(a,N) > 1:
            if pr: print('Factor found: %d'%gcd(a,N))
            if pr: print('COMPOSITE')
            return False
        r,s = 0,N-1
        while s % 2 == 0:
            r += 1
            s //= 2
        y = power_mod(a,s,N)
        if y == 1 or y == N-1: continue
        while r > 1:
            y = y*y % N
            if y == N-1: break
            r -= 1
        if r > 1: continue
        if pr: print('Miller-Rabin test failed with a = %d'%a)
        if pr: print('COMPOSITE')
        return False
    if pr: print('Miller-Rabin test passed with %d repeated tests.'%test_time)
    if pr: print('PRIME (probable)')
    return True


##### 3. Solovay-Strassen #####

def Solovay_Strassen(N:int,test_time=100,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    for _ in range(test_time):
        a = random.randint(2,N-1)
        if gcd(a,N) > 1:
            if pr: print('Factor found: %d'%gcd(a,N))
            if pr: print('COMPOSITE')
            return False
        if power_mod(a,N//2,N) != Jacobi(a,N) % N:
            if pr: print('Solovay-Strassen test failed with a = %d'%a)
            if pr: print('COMPOSITE')
            return False
    if pr: print('Solovay-Strassen test passed with %d repeated tests.'%test_time)
    if pr: print('PRIME (probable)')
    return True


##### 4. Lucas probable prime #####

def Lucas_prob_prime(N:int,P:int,Q:int,pr=True):
    assert N > 1, 'Tested integer should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    D = P*P - 4*Q
    if 1 < gcd(Q*D,N) < N:
        if pr: print('Factor found: %d'%gcd(Q*D,N))
        if pr: print('COMPOSITE')
        return False
    assert Q*D % N > 0, 'Tested integer should NOT DEVIDE Q*D'
    if Lucas_sequence(N-Jacobi(D,N),P,Q,N)[0] > 0:
        if pr: print('Lucas probable prime test failed.')
        if pr: print('COMPOSITE')
        return False
    if pr: print('Lucas probable prime test passed.')
    if pr: print('PRIME (probable)')
    return True

def Strong_Lucas_prob_prime(N:int,P:int,Q:int,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    D = P*P - 4*Q
    if 1 < gcd(Q*D,N) < N:
        if pr: print('Factor found: %d'%gcd(Q*D,N))
        if pr: print('COMPOSITE')
        return False
    assert Q*D % N > 0, 'Tested integer should NOT DEVIDE Q*D'
    r,s = 0,N-Jacobi(D,N)
    while s % 2 == 0:
        r += 1
        s //= 2
    U,V = Lucas_sequence(s,P,Q,N)
    if U*V == 0:
        if pr: print('Strong Lucas probable prime test passed with Us == 0 or Vs == 0 (mod N).')
        if pr: print('PRIME (probable)')
        return True
    while r > 1:
        U,V = U*V%N,D*U*U+V*V
        V = (V+V%2*N)//2 % N
        if V == 0:
            if pr: print('Strong Lucas probable prime test passed with V_(s*2^i) == 0 (mod N).')
            if pr: print('PRIME (probable)')
            return True
        r -= 1
    if pr: print('Strong Lucas probable prime test failed.')
    if pr: print('COMPOSITE')
    return False

def Extra_strong_Lucas_prob_prime(N:int,b:int,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    D = b*b - 4
    if 1 < gcd(D,N) < N:
        if pr: print('Factor found: %d'%gcd(D,N))
        if pr: print('COMPOSITE')
        return False
    assert D % N > 0, 'Tested integer should NOT DEVIDE D'
    r,s = 0,N-Jacobi(D,N)
    while s % 2 == 0:
        r += 1
        s //= 2
    U,V = Lucas_sequence(s,b,1,N)
    if U == 0 and (V == 2 or V == N-2) or V == 0:
        if pr: print('Extra strong Lucas probable prime test passed with (Us == 0 and Vs == ±2) or Vs == 0 (mod N).')
        if pr: print('PRIME (probable)')
        return True
    while r > 1:
        U,V = U*V%N,D*U*U+V*V
        V = (V+V%2*N)//2 % N
        if V == 0:
            if pr: print('Extra strong Lucas probable prime test passed with V_(s*2^i) == 0 (mod N).')
            if pr: print('PRIME (probable)')
            return True
        r -= 1
    if pr: print('Extra strong Lucas probable prime test failed.')
    if pr: print('COMPOSITE')
    return False


##### 5. Baillie-PSW #####

def Baillie_PSW(N:int,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    # trial division and perfert square test
    factor = 1
    for m in [2,3,5,7,11,13,17,19]:
        if m * m > N and factor == 1:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0: factor = m
    if powrt_int(N,2): factor = powrt_int(N,2)    # perfert square test
    if factor > 1:
        if pr: print('Factor found: %d'%powrt_int(N,2))
        if pr: print('COMPOSITE')
        return False
    ## Miller-Rabin based on 2
    r,s = 0,N-1
    while s % 2 == 0:
        r += 1
        s //= 2
    y = power_mod(2,s,N)
    if y != 1 and y != N-1:
        while r > 1:
            y = y*y % N
            if y == N-1: break
            r -= 1
        if r == 1:
            if pr: print('Miller-Rabin test (basis 2) failed.')
            if pr: print('COMPOSITE')
            return False
    if pr: print('Miller-Rabin test (basis 2) passed.')
    ## Find D
    D = 5
    while Jacobi(D,N) > -1: D = -(D+2) if D > 0 else -(D-2)
    ## strong Lucas probable prime test
    P,Q = 1,(1-D)//4
    if 1 < gcd(Q*D,N) < N:
        if pr: print('Factor found: %d'%gcd(Q*D,N))
        if pr: print('COMPOSITE')
        return False
    r,s = 0,N+1
    while s % 2 == 0:
        r += 1
        s //= 2
    U,V = Lucas_sequence(s,P,Q,N)
    SLPP = U*V == 0
    while r > 1:
        U,V = U*V%N,D*U*U+V*V
        V = (V+V%2*N)//2 % N
        if V == 0: SLPP = True
        r -= 1
    if not SLPP:
        if pr: print('Strong Lucas probable prime test failed with Us,Vs != 0 and V_(2^i*s) != 0 (mod N), and D,Q = %d,%d'%(D,Q))
        if pr: print('COMPOSITE')
        return False
    V = D*U*U+V*V
    V = (V+V%2*N)//2 % N
    if V == 2*Q%N and pr: print('Strong Lucas probable prime test passed with D,Q = %d,%d'%(D,Q))
    elif pr: print('Strong Lucas probable prime test failed with V_(N+1) != 2*Q (mod N), and D,Q = %d,%d'%(D,Q))
    if pr: print('PRIME (probable)' if V == 2*Q%N else 'COMPOSITE')
    return V == 2*Q%N


##### 6. Lucas #####

def Lucas(N:int,test_time=100,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    Qs = trial_div_factorize(N-1)
    if pr: print('Prime factors of N-1: %s'%Qs)
    for _ in range(test_time):
        a = random.randint(2,N-2)
        if gcd(a,N) > 1:
            if pr: print('Factor found: %d'%gcd(a,N))
            if pr: print('COMPOSITE')
            return False
        if power_mod(a,N-1,N) != 1:
            if pr: print('Fermat test failed with a = %d'%a)
            if pr: print('COMPOSITE')
            return False
        q_test = True
        for q in Qs:
            if power_mod(a,(N-1)//q,N) == 1:
                q_test = False
                break
        if q_test:
            if pr: print('Lucas test passed with a = %d'%a)
            if pr: print('PRIME')
            return True
    if pr: print('Lucas test failed with %d repeated tests.'%test_time)
    if pr: print('COMPOSITE (probable)')
    return False


##### 7. Pocklington #####

def Pocklington(N:int,test_time=100,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    for m in [2,3,5,7,11,13,17,19]:
        if m * m > N:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0:
            if pr: print('Factor found: %d'%m)
            if pr: print('COMPOSITE')
            return False
    F = FR_factorize(N-1,balance=True)[0]
    if pr: print('Factorization of N-1: F = %d, R = %d'%(F,(N-1)//F))
    Qs = trial_div_factorize(F)
    if pr: print('Prime factors of F: %s'%Qs)
    for q in Qs:
        found_a = False
        for _ in range(test_time):
            a = random.randint(2,N-2)
            if gcd(a,N) > 1:
                if pr: print('Factor found: %d'%gcd(a,N))
                if pr: print('COMPOSITE')
                return False
            if power_mod(a,N-1,N) != 1:
                if pr: print('Fermat test failed with a = %d'%a)
                if pr: print('COMPOSITE')
                return False
            x = power_mod(a,(N-1)//q,N) - 1
            if gcd(x,N) == 1:
                found_a = True
                break
        if not found_a:
            if pr: print('Pocklington test failed with no basis found for q = %d'%q)
            if pr: print('COMPOSITE (probable)')
            return False
    if pr: print('Pocklington test passed.')
    if pr: print('PRIME')
    return True


##### 8. Proth (for P=A*2^n+1) & Pépin (for F_n=2^{2^n}+1) #####

def Proth(A:int,n:int,test_time=100,pr=True):
    assert n > 0, 'Integer n should be POSITIVE'
    assert A > 0 and A % 2 > 0 and A < 2**n, 'Integer A should be ODD and LESS than 2^n.'
    P = A * 2**n + 1
    if pr: print('P = %d'%P)
    for _ in range(test_time):
        a = random.randint(2,P-1)
        if gcd(a,P) > 1:
            if pr: print('Factor found: %d'%gcd(a,P))
            if pr: print('COMPOSITE')
            return False
        y = power_mod(a,A,P)
        for _ in range(n-1): y = y*y % P
        if y == P-1:
            if pr: print('Proth test passed with a = %d'%a)
            if pr: print('PRIME')
            return True
        if y == 1 == -Jacobi(a,P) or y != 1:
            if pr: print('Fermat or Solovay-Strassen test failed with a = %d'%a)
            if pr: print('COMPOSITE')
            return False
    if pr: print('Proth test failed with %d repeated tests.'%test_time)
    if pr: print('COMPOSITE (probable)')
    return False

def Pepin(n,pr=True):
    assert n >= 0, 'Integer n should be NOT NEGATIVE.'
    m = 2**n
    Fn = 2**m + 1
    if pr: print('Fn = %d'%Fn)
    if n == 0:
        if pr: print('PRIME')
        return True
    y = 3
    for _ in range(m-1): y = y*y % Fn
    if pr: print('3^((Fn-1)/2) %s= -1 (mod Fn)'%('=' if y==Fn-1 else '!'))
    if pr: print('PRIME' if y==Fn-1 else 'COMPOSITE')
    return y == Fn-1


##### 9. Factorizing Lucas probable prime #####

def Factorizing_Lucas_prob_prime_v1(N:int,test_time=100,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    factor = 1
    for m in [2,3,5,7,11,13,17,19]:         # trial division
        if m * m > N and factor == 1:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0: factor = m
    if powrt_int(N,2): factor = powrt_int(N,2)
    if factor > 1:
        if pr: print('Factor found: %d'%factor)
        if pr: print('COMPOSITE')
        return False
    Qs = trial_div_factorize(N+1)
    if pr: print('Prime factors of N+1: %s'%Qs)
    D = 5
    while test_time:
        P,Q = 1,(1-D)//4
        if Jacobi(D,N) > -1 or Q*D % N == 0:
            D = -(D+2) if D > 0 else -(D-2)
            continue
        if 1 < gcd(Q*D,N) < N:
            if pr: print('Factor found: %d'%gcd(Q*D,N))
            if pr: print('COMPOSITE')
            return False
        U,V = Lucas_sequence(N+1,P,Q,N)
        if U > 0 or V != 2*Q%N:
            if pr: print('Lucas probable prime test failed with D,Q = %d,%d'%(D,Q))
            if pr: print('COMPOSITE')
            return False
        q_test = True
        for q in Qs:
            if Lucas_sequence((N+1)//q,P,Q,N)[0] == 0:
                q_test = False
                break
        if q_test:
            if pr: print('Factorizing Lucas probable prime test passed with D,Q = %d,%d'%(D,Q))
            if pr: print('PRIME')
            return True
        D = -(D+2) if D > 0 else -(D-2)
        test_time -= 1
    if pr: print('Factorizing Lucas probable prime test failed with %d repeated tests.'%test_time)
    if pr: print('COMPOSITE (probable)')
    return False

def Factorizing_Lucas_prob_prime_v2(N:int,test_time=100,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    factor = 1
    for m in [2,3,5,7,11,13,17,19]:
        if m * m > N and factor == 1:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0: factor = m
    if powrt_int(N,2): factor = powrt_int(N,2)
    if factor > 1:
        if pr: print('Factor found: %d'%factor)
        if pr: print('COMPOSITE')
        return False
    for _ in range(10):
        F = FR_factorize(N+1,balance=True)[0]
        if (F-1)**2 > N: break
    if (F-1)**2 <= N: F = N
    if pr: print('Factorization of N+1: F = %d, R = %d'%(F,(N+1)//F))
    Qs = trial_div_factorize(F)
    if pr: print('Prime factors of F: %s'%Qs)
    for q in Qs:
        found_PQ = False
        D = 5
        while test_time:
            P,Q = 1,(1-D)//4
            if Jacobi(D,N) > -1 or Q*D % N == 0:
                D = -(D+2) if D > 0 else -(D-2)
                continue
            if 1 < gcd(Q*D,N) < N:
                if pr: print('Factor found: %d'%gcd(Q*D,N))
                if pr: print('COMPOSITE')
                return False
            U,V = Lucas_sequence(N+1,P,Q,N)
            if U > 0 or V != 2*Q%N:
                if pr: print('Lucas probable prime test failed with D,Q = %d,%d'%(D,Q))
                if pr: print('COMPOSITE')
                return False
            u = Lucas_sequence((N+1)//q,P,Q,N)[0]
            if gcd(u,N) == 1:
                found_PQ = True
                break
            D = -(D+2) if D > 0 else -(D-2)
            test_time -= 1
        if not found_PQ:
            if pr: print('Factorizing Lucas probable prime test failed with no P,Q found for q = %d'%q)
            if pr: print('COMPOSITE (probable)')
            return False
    if pr: print('Factorizing Lucas probable prime test passed.')
    if pr: print('PRIME')
    return True


##### 10. Lucas-Lehmer (for M_p=2^p-1) & Lucas-Lehmer-Riesel (for N=A*2^n-1) #####

def Lucas_Lehmer(p:int,pr=True):
    assert p > 1, 'Prime p should be NO LESS THAN 2'
    Mp = 2**p - 1
    if pr: print('Mp = %d'%Mp)
    if p == 2:
        if pr: print('PRIME')
        return True
    # trial test for p
    factor = 1
    for m in [2,3,5]:
        if m * m > p: break
        if p % m == 0: factor = m
    if factor == 1 and m * m < p:
        A,i,m = 0,1,7
        B = [1,7,11,13,17,19,23,29]
        while m * m <= p:
            if p % m == 0:
                factor = m
                break
            i = (i+1) % 8
            if i == 0: A += 30
            m = A + B[i]
    if factor > 1:
        if pr: print('Factor found: %d'%(2**factor-1))
        if pr: print('COMPOSITE')
        return False
    # Lucas_Lehmer for Mp
    S = 4 % Mp
    for _ in range(p-2): S = (S*S-2) % Mp
    if pr: print('Sp %s= 0 (mod Mp)'%('!' if S > 0 else '='))
    if pr: print('COMPOSITE' if S > 0 else 'PRIME')
    return S == 0

def Lucas_Lehmer_Riesel(A:int,n:int,pr=True):
    assert A % 2 > 0 and 0 < A < 2**n, 'A should be POSITIVE ODD less than 2^n'
    assert n >= 0, 'n should be NOT NEGTIVE'
    N = A * 2**n - 1
    if pr: print('N = %d'%N)
    if A % 3 > 0: P = 4
    else:
        P = 5
        while Jacobi(P-2,N) != 1 or Jacobi(P+2,N) != -1:
            if P == 5: P = 8
            elif P == 9: P = 11
            else: P += 1
    S = Lucas_sequence(A,P,1,N)[1]
    for _ in range(n-2): S = (S*S-2) % N
    if pr: print('Sn %s= 0 (mod N)'%('!' if S > 0 else '='))
    if pr: print('COMPOSITE' if S > 0 else 'PRIME')
    return S == 0


##### 11. ECPP: Goldwasser-Kilian #####

def Goldwasser_Kilian(N:int,test_time=50,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    test_number = 1
    while test_number <= test_time:
        if pr: print('Test No.%d'%test_number)
        ## Step 1: randomly choose a,b,x,y
        a,x,y = [random.randint(0,N-1) for _ in range(3)]
        b = (y*y - x**3 - a*x) % N
        if (4*a**3 + 27*b*b) % N == 0: continue
        P = (x,y)
        if pr: print('a,b = %d,%d   P = %s'%(a,b,P))
        test_number += 1
        ## Step 2: Schoof's algorithm to get order m
        ep = Elliptic_curve(N,a,b)
        if pr: print('running Schoof\'s algorithm ...')
        m = ep.schoof()
        if type(m) == dict:
            if pr: print('\n%s'%m['info'])
            if pr: print('COMPOSITE')
            return False
        if pr: print('Order of group: m = %d'%m)
        ## Step 3: Find factor s of m, and factorize s
        for _ in range(10):
            s = FR_factorize(m,balance=True)[0]
            if s > (N**(1/4)+1)**2: break
        if s <= (N**(1/4)+1)**2: s = m
        if pr: print('factor of m: s = %d'%s)
        Q1,Q2 = prob_prime_factorize(s)
        if pr: print('factors of s: %s, %s'%(Q1,Q2))
        ## Steo 4,5,6: Compute mP and (m/q)P
        Pm = ep.times(m,P)
        if type(Pm) == dict:
            if pr: print('\n%s'%Pm['info'])
            if pr: print('COMPOSITE')
            return False
        if Pm != (0,):
            if pr: print('\nGroup order test failed: mP != 0')
            if pr: print('COMPOSITE')
            return False
        if pr: print('mP = 0')
        q_test = True
        for q in Q1|Q2:
            Pm = ep.times(m//q,P)
            if type(Pm) == dict:
                if pr: print('\n%s'%Pm['info'])
                if pr: print('COMPOSITE')
                return False
            if pr: print('(m/%d)P = %s'%(q,Pm))
            if Pm == (0,):
                if pr: print('Let\'s try again ...\n')
                q_test = False
                break
        if not q_test: continue
        ## Step 7: Recursion primality test for all q in Q2
        q_test = True
        for q in Q2:
            if pr: print('\nRecursion: ECPP for q = %d ...'%q)
            if not Goldwasser_Kilian(q,pr=pr):
                if pr: print('\nq = %d is composite, try again ...\n'%q)
                q_test = False
                break
            elif pr: print('\nq = %d is prime.'%q)
        if not q_test: continue
        if pr: print('\nECPP success with a,b,x,y = %d,%d,%d,%d'%(a,b,x,y))
        if pr: print('PRIME')
        return True
    ## Step 8: Sufficient test to tell N is composite
    if pr: print('\nElliptic curve with proper order not found.')
    if pr: print('COMPOSITE (probable)')
    return False


##### 12. ECPP: Atkin-Morain #####

def Atkin_Morain(N:int,max_D=500,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    square_factors = [i*i for i in range(2,int(math.sqrt(max_D))+1)]
    test_number = 1
    for D in random.sample(range(1,max_D+1),max_D):
        ## Step 1: choose D, and solve equation 4N == u^2 + Dv^2
        D = random.randint(1,max_D)
        if 1 < gcd(D,N) < N:
            if pr: print('Test No.%d:  D = %d\n'%(test_number,D))
            if pr: print('Factor found: %d'%gcd(D,N))
            if pr: print('COMPOSITE')
            return False
        if any(D % s == 0 for s in square_factors): continue
        if Jacobi(-D,N) < 1: continue
        if pr: print('Test No.%d:  D = %d'%(test_number,D))
        test_number += 1
        RCP = reduced_class_polynomial(D)
        if type(RCP) == dict:
            if pr: print(RCP['info'],'\n')
            continue
        RCP.modulo = N
        uvs = Cornacchia_4N(N,D)
        if type(uvs) == dict:
            if pr: print('\n%s'%uvs['info'])
            if pr: print('COMPOSITE')
            return False
        if not uvs:
            if pr: print('Solution not found for u^2 + Dv^2 = 4N ...\n')
            continue
        # list all posible order m
        if len(uvs) == 1:
            u = list(uvs)[0][0]
            ms = [N+1+u,N+1-u]
        else:
            ms = []
            for u,_ in uvs:
                if u % 2 == 0: ms = ms + [N+1+u,N+1-u]
        if pr: print('u = %d, probable values of order m: %s'%(u,ms))
        ## Step 2: complex multiplication method to get a,b, and randomly choose a point P on the curve
        a,b = complex_multiplication(N,D,RCP,u%2)
        if pr: print('a = %d, b = %d'%(a,b))
        ep = Elliptic_curve(N,a,b)
        P = ep.random_point()
        if type(P) == dict:
            if pr: print('\n%s'%P['info'])
            if pr: print('COMPOSITE')
            return False
        if pr: print('randomly choose a point: P(x,y) = %s'%(P,))
        ## Step 3: find m such that mP == 0
        for m in ms+[0]:
            Pm = ep.times(m,P)
            if type(Pm) == dict:                        # information received which tells N is composite
                if pr: print('\n%s'%Pm['info'])
                if pr: print('COMPOSITE')
                return False
            if Pm == (0,): break                        # mP == 0
        if m == 0:                                               # all mP != 0, return composite
            if pr: print('\nGroup order test failed: mP != 0')
            if pr: print('COMPOSITE')
            return False
        if pr: print('Order of group: m = %d'%m)
        ## Step 4: Find factor s of m, and factorize s
        for _ in range(10):
            s = FR_factorize(m,balance=True)[0]
            if s > (N**(1/4)+1)**2: break
        if s <= (N**(1/4)+1)**2: s = m
        if pr: print('factor of m: s = %d'%s)
        Q1,Q2 = prob_prime_factorize(s)
        if pr: print('factors of s: %s, %s'%(Q1,Q2))
        ## Steo 5,6: Compute (m/q)P
        q_test = True
        for q in Q1|Q2:
            Pm = ep.times(m//q,P)
            if type(Pm) == dict:                    # information received which tells N is composite
                if pr: print('\n%s'%Pm['info'])
                if pr: print('COMPOSITE')
                return False
            if pr: print('(m/%d)P = %s'%(q,Pm))
            if Pm == (0,):                          # (m/q)P == 0, back to step 1 and choose another D
                if pr: print('Let\'s try again ...\n')
                q_test = False
                break
        if not q_test: continue
        ## Step 7: Recursion primality test for all q in Q2
        q_test = True
        for q in Q2:
            if pr: print('\nRecursion: ECPP for q = %d ...'%q)
            if not Atkin_Morain(q,pr=pr):          # one of the probable prime factors is said to be composite
                if pr: print('\nq = %d is composite.'%q)
                q_test = False
                break
            elif pr: print('\nq = %d is prime, try again ...\n'%q)
        if not q_test: continue
        if pr: print('\nECPP success with m,a,b,x,y = %d,%d,%d,%d,%d'%(m,a,b,P[0],P[1]))    # All probable prime factors are definitely prime, so N is prime
        if pr: print('PRIME')
        return True
    ## Step 8: Sufficient test to tell N is composite
    if pr: print('\nElliptic curve with proper order not found.')
    if pr: print('COMPOSITE (probable)')
    return False


##### 13. ECPP: Tsumura #####

def Tsumura(k:int,n:int,pr=True):
    assert k >= 2, 'k should be no less than 2.'
    assert n >= 1 and n&1, 'n should be POSITIVE and ODD.'
    N = 2**k * n - 1
    M = (N**0.25+1)**2
    assert n < N/M or n > M, 'n not in applicable range.'
    if pr: print('N = %d'%N)
    while True:
        x = random.randint(1,N-1)
        jacobi = Jacobi(x,N)
        if jacobi == -1: break
        if jacobi == 0:
            if pr: print('Factor found: %d'%gcd(x,N))
            if pr: print('COMPOSITE')
            return False
    if pr: print('randomly choose x = %d'%x)
    while True:
        while True:
            y = random.randint(1,N-1)
            m = (x**3-y*y) % N
            if m == 0: continue
            jacobi = Jacobi(m,N)
            if jacobi == 1: break
            if jacobi == 0:
                if pr: print('Factor found: %d'%gcd(m,N))
                if pr: print('COMPOSITE')
                return False
        m = m * exgcd(x,N)[1] % N
        ep = Elliptic_curve(N,-m,0)
        P = (x,y)
        if pr: print('randomly choose y = %d\nm = %d'%(y,m))
        if n < N/M:
            if pr: print('Range of n:  n < N/(N^(1/4)+1)^2')
            R = ep.times(n,P)
            if type(R) == dict:
                if pr: print(R['info'])
                if pr: print('COMPOSITE')
                return False
            if R == (0,):
                if pr: print('Order test failed: nP == 0')
                if pr: print('COMPOSITE')
                return False
            X = R[0]
            for i in range(1,k):
                S = X**3 - m*X
                if S % N == 0:
                    if pr: print('Order test failed: (2^%d)nP == 0'%i)
                    if pr: print('COMPOSITE')
                    return False
                if gcd(S,N) > 1:
                    if pr: print('Factor found: %d'%gcd(S,N))
                    if pr: print('COMPOSITE')
                    return False
                X = (X*X+m)**2 * exgcd(4*S,N)[1] % N
            S = (X**3 - m*X) % N
            if pr: print('Order test failed: (2^%d)nP != 0'%k if S > 0 else 'ECPP success with m,x,y = %d,%d,%d'%(m,x,y))
            if pr: print('COMPOSITE' if S > 0 else 'PRIME')
            return not S
        if pr: print('Range of n:  n > (N^(1/4)+1)^2')
        R = P
        for _ in range(k):
            R = ep.plus(R,R)
            if type(R) == dict:
                if pr: print(R['info'])
                if pr: print('COMPOSITE')
                return False
        if pr: print('(2^k)P = %s'%(R,))
        if R == (0,):
            if pr: print('try another y ...')
            continue
        Qs = trial_div_factorize(n)
        if pr: print('prime factors of n: %s'%Qs)
        q_test = True
        for q in Qs:
            Rn = ep.times(n//q,R)
            if type(Rn) == dict:                    # information received which tells N is composite
                if pr: print(Rn['info'])
                if pr: print('COMPOSITE')
                return False
            if pr: print('(2^k)(n/%d)P = %s'%(q,Rn))
            if Rn == (0,):
                if pr: print('try another y ...')
                q_test = False
                break
        if not q_test: continue
        Rn = ep.times(n,R)
        if pr: print('ECPP success with m,x,y = %d,%d,%d'%(m,x,y) if Rn == (0,) else 'Order test failed: (2^k)nP != 0')
        if pr: print('PRIME' if Rn == (0,) else 'COMPOSITE')
        return Rn == (0,)


##### 14. High-order field & cyclotomic polynomial #####

def Quartic_field(N:int,max_D=120,max_C=120,test_time=100,pr=True):
    assert N > 1, 'Tested integer should be NO LESS THAN 2'
    factor = 1
    for m in [2,3,5,7,11,13,17,19]:                     # trial division
        if m * m > N and factor == 1:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0: factor = m
    if powrt_int(N,2) > 0: factor = powrt_int(N,2)      # perfect square test
    if factor > 1:
        if pr: print('Factor found: %d'%factor)
        if pr: print('COMPOSITE')
        return False
    if pr: print('Phi = N*N + 1 = %d'%(N*N+1))
    F4,R4 = 1,N*N//2+1
    while (F4+1)**3 <= N*N:
        B4 = (N**(2/3)+1) / F4
        s = Pollard_rho(R4,test_time=10,loop_time=max(int(math.sqrt(B4)),2**16))
        if s == 1: break
        R4 //= s
        r = gcd(s,R4)
        while r > 1:
            s *= r
            R4 //= r
            r = gcd(s,R4)
        F4 *= s
    if pr: print('F4 = %d, R4 = %d'%(F4,R4))
    Q4 = trial_div_factorize(F4)
    if R4 > 1: Q4.add(R4)
    if pr: print('All factors to be tested: %s'%Q4)
    D = 5
    while abs(D) < max_D:
        if 1 < gcd(D,N) < N:
            if pr: print('Factor found: %d'%gcd(D,N))
            if pr: print('COMPOSITE')
            return False
        if D % N == 0 or Jacobi(D,N) == 1:
            D = -(D+2) if D > 0 else -(D-2)
            continue
        for C in range(max_C):
            if 1 < gcd(C*C-16*D,N) < N:
                if pr: print('Factor found: %d'%gcd(C*C-16*D,N))
                if pr: print('COMPOSITE')
                return False
            if (C*C-16*D) % N == 0 or Jacobi(C*C-16*D,N) == 1: continue
            if pr: print('\nSet parameter: D,C = %d,%d'%(D,C))
            q_test = True
            for q in Q4:
                param_found = False
                for _ in range(test_time):
                    H,K = [random.randint(0,N-1) for _ in [0,1]]
                    if 1 < gcd(H*H-K*K*D,N) < N:
                        if pr: print('Factor found: %d'%gcd(H*H-K*K*D,N))
                        if pr: print('COMPOSITE')
                        return False
                    if (H*H-K*K*D) % N == 0: continue
                    p = (2*H*H + H*K*C + 2*K*K*D) % N
                    P1,P2,Q = 4*p % N, 4*(p*p-D) % N, (p*p-H*H*C-8*H*K*D-K*K*C*D+D) % N
                    if V2_sequence(N*N+1,P1,P2,Q,N)[1] % N > 0:
                        if pr: print('Quartic field test failed: V_(N^2+1) != 0 (mod N), with D,C,H,K = %d,%d,%d,%d'%(D,C,H,K))
                        if pr: print('COMPOSITE')
                        return False
                    V = V2_sequence((N*N+1)//q,P1,P2,Q,N)[1]
                    if gcd(V,N) == 1:
                        param_found = True
                        break
                if not param_found:
                    q_test = False
                    break
            if q_test:
                if pr: print('Quartic field test passed with D,C = %d,%d'%(D,C))
                if pr: print('PRIME')
                return True
        D = -(D+2) if D > 0 else -(D-2)
    if pr: print('Proper parameters not found.')
    if pr: print('COMPOSITE (probable)')
    return False

def Cubic_Sextic_field(N:int,theta:int,max_G=120,test_time=100,pr=True):
    assert N > 1, 'Tested integer should be NO LESS THAN 2'
    assert abs(theta) == 1, 'Parameter theta must be 1 or -1.'
    factor = 1
    for m in [2,3,5,7,11,13,17,19]:                     # trial division
        if m * m > N and factor == 1:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0: factor = m
    if powrt_int(N,2) > 0: factor = powrt_int(N,2)      # perfect square test
    if factor > 1:
        if pr: print('Factor found: %d'%factor)
        if pr: print('COMPOSITE')
        return False
    Phi = N*N + theta*N + 1
    if pr: print('Phi = N^2 %s N + 1 = %d'%('+' if theta>0 else '-',Phi))
    F,R = 1,Phi
    while F**3 <= N*N:
        s = Pollard_rho(R,test_time=10,loop_time=max(int(R**(1/4)),1<<16))
        if s == 1: s = R
        R //= s
        r = gcd(s,R)
        while r > 1:
            s *= r
            R //= r
            r = gcd(s,R)
        F *= s
    if pr: print('F = %d, R = %d'%(F,R))
    Qs = trial_div_factorize(F)
    if pr: print('All factors to be tested: %s'%Qs)
    for P,S,T in [(7,1,1),(19,7,1),(31,4,2),(61,1,3),(73,7,3),(97,19,1),(103,13,3),(151,19,3),(163,25,1),
                  (181,7,5),(211,13,5),(223,28,2),(229,22,4),(307,16,6),(331,1,7),(349,37,1),(373,13,7),
                  (397,34,4),(409,31,5),(421,19,7),(439,28,6),(457,10,8),(487,25,7)]:
        if 1 < gcd(P,N) < N:
            if pr: print('Factor found: %d'%gcd(P,N))
            if pr: print('COMPOSITE')
            return False
        if N == P:
            if pr: print('Contained in prime list.')
            if pr: print('PRIME')
            return True
        if power_mod(N%P,(P-1)//3,P) == 1: continue
        for G in range(1,max_G):
            if 1 < gcd(G,N) < N:
                if pr: print('Factor found: %d'%gcd(G,N))
                if pr: print('COMPOSITE')
                return False
            if G % N == 0 or Jacobi(G,N) == -theta: continue
            if pr: print('\nSet parameter: P,S,G = %d,%d,%d'%(P,S,G))
            q_test = True
            for q in Qs:
                param_found = False
                for _ in range(test_time):
                    H,K,L = [random.randint(0,N-1) for _ in [0,1,2]]
                    a,b = 3*P,P*S
                    m1,m2,m3 = (H*H+2*b*L*K-G) % N, (b*L*L+2*H*K+2*a*K*L) % N, (K*K+a*L*L+2*H*L) % N
                    d1,d2,d3 = ((m1+a*m3)**2 - m2*(a*m2+b*m3)) % N, (b*m3*m3-m1*m2) % N, (m2*m2-m1*m3-a*m3*m3) % N
                    r = (d1*m1 + b*d2*m3 + b*d3*m2) % N
                    A,B,C = (4*G*d1+2*r) % N, 4*G*d2 % N, 4*G*d3 % N
                    P1,P2,P3,Q = (3*A+2*a*C) % N, ((A+a*C)*(3*A+a*C)-a*B*B-3*b*B*C) % N, (A*(A+a*C)**2+b*B**3+b*b*C**3-3*b*A*B*C-a*A*B*B-a*b*B*C*C) % N, r*r % N
                    Delta = (P1*P2*(P1*P2+18*P3) - 4*P2**3 - P3*(4*P1**3+27*P3))%N
                    E = ((P3+4*Q*P1)**2 - 4*Q*(P2+4*Q)**2) % N
                    if 1 < gcd(Delta*E*Q,N) < N:
                        if pr: print('Factor found: %d'%gcd(Delta*E*Q,N))
                        if pr: print('COMPOSITE')
                        return False
                    if (Delta*E*Q) % N == 0: continue
                    _,v1,v2 = V3_sequence(Phi,P1,P2,P3,Q,N)
                    if gcd(v1,v2) % N > 0:
                        if pr: print('%s field test failed: C_(N^2%sN+1) != 0 mod N, with P,S,G,H,K,L = %d,%d,%d,%d,%d,%d'%((('Cubic','+') if theta>0 else ('Sextic','-'))+(P,S,G,H,K,L)))
                        if pr: print('COMPOSITE')
                        return False
                    _,v1,v2 = V3_sequence(Phi//q,P1,P2,P3,Q,N)
                    if gcd(v1,v2,N) == 1:
                        param_found = True
                        break
                if not param_found:
                    q_test = False
                    break
            if q_test:
                if pr: print('%s field test success with P,S,G = %d,%d,%d'%('Cubic' if theta>0 else 'Sextic',P,S,G))
                if pr: print('PRIME')
                return True
    if pr: print('Parameters not found.')
    if pr: print('COMPOSITE (probable)')
    return False


##### 15. APR-CL test #####

def APR_CL(N:int,pr=True):
    assert N > 1, 'Tested integer should be NO LESS THAN 2'
    for m in [2,3,5,7,11,13,17,19]:     # trial division
        if m * m > N:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0:
            if pr: print('Factor found: %d'%m)
            if pr: print('COMPOSITE')
            return False
    ## Step 1
    for t in [2,6,28,20,12,24,40,48,36,60,72,120,300,180,360,420,720,840,1080,1200,
              1260,1680,2520,7920,5040,12600,15120,25200,55440,113400,110880,221760,504000,720720,
              1441440,4324320,24504480,73513440,367567200,1396755360,6983776800]:      # A list of parameter t, which supports tests for at most 6011 digits
        s = e(t)
        if s * s > N: break
    if pr: print('Step 1\nt,s = %d,%d'%(t,s))
    if gcd(t*s,N) > 1:
        if pr: print('\nFactor found: %d'%gcd(s*t,N))
        if pr: print('COMPOSITE')
        return False
    Euclidean = trial_div_factorize(s)                     # list all Euclidean primes q
    if pr: print('Euclidean primes: %s'%Euclidean)
    ## Step 2
    if pr: print('\nStep 2')
    js,js_ast,js_tag = {},{},{}
    for q in Euclidean:
        if q == 2: continue
        f = generate_f(q)                    # Step 2.1
        for p in trial_div_factorize(q-1):   # Step 2.2
            k = v(p,q-1)                        # Step 2.2.1
            if p**k > 2: js[(p,q)] = calculate_jpq(p,k,q,f)    # Step 2.2.2
            if p == 2 and k > 2: js_ast[q],js_tag[q] = calculate_j2q(k,q,js[(2,q)],f)   # Step 2.2.3
        if pr: print('j(p,q) for q = %d are calculated.'%q)
    ## Stage B: loop of Step 3,4,5
    if pr: print('\nStep 3,4,5')
    for p in trial_div_factorize(t):
        lbd = p > 2 and power_mod(N,p-1,p*p) != 1       # Step 3
        for q in Euclidean:                           # Step 4
            if (q-1) % p > 0: continue
            k = v(p,q-1)                              # Step 4.1
            uk,vk = N//p**k,N%p**k
            if p > 2: J0,Jv = calculate_Jpq(p,k,vk,js[(p,q)],N)                 # Step 4.2
            elif k == 1: J0,Jv = q,1
            elif k == 2: J0,Jv = js[(p,q)]**2*q, 1 if vk==1 else js[(p,q)]**2
            else: J0,Jv = calculate_J2q(k,vk,js_ast[q],js_tag[q],N)
            z = power_mod(J0,uk,N) * Jv               # Step 4.3
            h = find_h(p,k,z,N)                       # find h
            if h == p**k:
                if pr: print('\nAPR-CL test failed at Step 4.3 with p,q = %d,%d'%(p,q))
                if pr: print('COMPOSITE')
                return False
            if h % p > 0 and (p > 2 or k == 1 and N % 4 == 1): lbd = True            # Step 4.4
            if h % 2 > 0 and p == 2 and k > 1 and not lbd:                    # Step 4.5
                if power_mod(q,N//2,N) < N-1:
                    if pr: print('\nAPR-CL test failed at Step 4.5 with p,q = %d,%d'%(p,q))
                    if pr: print('COMPOSITE')
                    return False
                lbd = True
        if lbd and pr: print('lambda_(p=%d) = True'%p)
        if not lbd:                                           # Step 5
            if pr: print('lambda_(p=%d) = False ...'%p)
            q = 2*p+1 if p>2 else 5                        # Step 5.1
            try_time = 0
            while try_time < len(Euclidean)+100:              # max attempt times
                if trial_division_30(q):                 # simply trial division
                    if q not in Euclidean and power_mod(N,(q-1)//p,q) != 1: break
                    try_time += 1
                q += 2*p if p>2 else (4 if q%6==1 else (8 if N%4==3 else 2))
            if try_time == len(Euclidean)+100:               # q' not found within max attempt times, algorithm fails
                if pr: print('\nUnable to test primality as q\' not found.')
                return
            if pr: print('prime q\' = %d is found'%q)
            if N % q == 0:                  # Step 5.2
                if pr: print('\nFactor found: %d'%q)
                if pr: print('COMPOSITE')
                return False
            k = 2 if p==2 and N%4==3 else 1       # Step 5.3
            uk,vk = N//p**k,N%p**k
            f = generate_f(q)
            if p**k > 2: j = calculate_jpq(p,k,q,f)
            if p > 2: J0,Jv = calculate_Jpq(p,1,vk,j,N)
            elif k == 1: J0,Jv = q,1
            else: J0,Jv = j*j*q, j*j
            if pr: print('j(p,q\') for q\' = %d are calculated.'%q)
            z = power_mod(J0,uk,N) * Jv      # Step 5.4
            h = find_h(p,k,z,N)
            if h % p == 0:
                if pr: print('\nAPR-CL test failed at Step 5.4 with p,q\' = %d,%d'%(p,q))
                if pr: print('COMPOSITE')
                return False
            if k == 2 and power_mod(q,N//2,N) < N-1:       # Step 5.5
                if pr: print('\nAPR-CL test failed at Step 5.5 with p,q\' = %d,%d'%(p,q))
                if pr: print('COMPOSITE')
                return False
            if pr: print('lambda_(p=%d) -> True'%p)
    ## Stage C: Step 6,7
    if pr: print('\nStep 6,7')
    r = 1
    for _ in range(1,t+1):
        r = r*N % s
        if r == 1:
            if pr: print('\nAPR-CL test success.')
            if pr: print('PRIME')
            return True
        if N % r == 0:
            if pr: print('\nFactor found: %d'%r)
            if pr: print('COMPOSITE')
            return False


##### 16. Frobenius probable prime test #####

def Frobenius_prob_prime(N:int,f:Polynomial,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    f.modulo = N
    d = f.degree
    Delta = f.discriminant()
    assert f.coef[0]*Delta % N > 0, 'N should not devide f(0)*Delta'
    if 1 < gcd(f.coef[0]*Delta,N) < N:
        if pr: print('Factor found: %d'%gcd(f.coef[0]*Delta,N))
        if pr: print('COMPOSITE')
        return False
    x = Polynomial(modulo=N)
    fk = f.copy()
    Nk = 1
    S = 0
    for k in range(1,d+1):
        Nk *= N
        Fk = poly_pow_mod(x,Nk,fk) - x
        if type(Fk) == Polynomial: Fk = gcd_poly(Fk,fk)
        if type(Fk) == Polynomial: Fk = Fk.monic()
        if type(Fk) == dict:                               # information received which tells N is composite
            if pr: print(Fk['info'])
            if pr: print('COMPOSITE')
            return False
        y = poly_pow_mod(x,N,Fk)
        if type(y) == Polynomial: y = poly_nest(Fk,y) % Fk
        if type(y) == dict:                                    # information received which tells N is composite
            if pr: print(y['info'])
            if pr: print('COMPOSITE')
            return False
        if y.coef:
            if pr: print('Frobenius step failed with k = %d'%k)
            if pr: print('COMPOSITE')
            return False
        fk //= Fk
        if k % 2 == 0: S += Fk.degree // k
    if fk.degree > 1:
        if pr: print('Factorization Step failed: fd(x) != 1')
        if pr: print('COMPOSITE')
        return False
    if (S%2, Jacobi(Delta,N)) not in [(0,1),(1,-1)]:
        if pr: print('Jacobi step failed.')
        if pr: print('COMPOSITE')
        return False
    if pr: print('Frobenius probable prime test passed.')
    if pr: print('PRIME (probable)')
    return True

def Strong_Frobenius_prob_prime(N:int,f:Polynomial,pr=True):         # An additional piece of code attached to function 'Frobenius_prob_prime'
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    if N <= 3:
        if pr: print('PRIME')
        return True
    if N % 2 == 0:
        if pr: print('Factor found: 2')
        if pr: print('COMPOSITE')
        return False
    f.modulo = N
    d = f.degree
    Delta = f.discriminant()
    assert f.coef[0]*Delta % N, 'N should not devide f(0)*Delta'
    if 1 < gcd(f.coef[0]*Delta,N) < N:
        if pr: print('Factor found: %d'%gcd(f.coef[0]*Delta,N))
        if pr: print('COMPOSITE')
        return False
    x = Polynomial(modulo=N)
    fk = f.copy()
    Nk = 1
    S = 0
    for k in range(1,d+1):
        Nk *= N
        Fk = poly_pow_mod(x,Nk,fk) - x
        if type(Fk) == Polynomial: Fk = gcd_poly(Fk,fk)
        if type(Fk) == Polynomial: Fk = Fk.monic()
        if type(Fk) == dict:
            if pr: print(Fk['info'])
            if pr: print('COMPOSITE')
            return False
        y = poly_pow_mod(x,N,Fk)
        if type(y) == Polynomial: y = poly_nest(Fk,y) % Fk
        if type(y) == dict:
            if pr: print(y['info'])
            if pr: print('COMPOSITE')
            return False
        if y.coef:
            if pr: print('Frobenius step failed with k = %d'%k)
            if pr: print('COMPOSITE')
            return False
        ## The additional code begins here
        r,s = 1,Nk-1
        while s % 2 == 0:
            r += 1
            s //= 2
        y = poly_pow_mod(x,s,Fk)
        if type(y) == Polynomial:
            yi = gcd_poly(y-1,Fk)               # 计算 Fk0(x)
            if type(yi) == Polynomial: yi = yi.monic()   # 首项归一化
            Fki = [yi]
            if type(yi) == dict: info = yi['info']
            else: info = 'Sqrt step failed with deg(Fki) mod k != 0' if yi.degree % k > 0 else ''
            # 'info' contains informations that N is composite. 'info' is empty if N is not found to be composite
            while len(Fki) <= r and not info:
                yi = gcd_poly(y+1,Fk)
                if type(yi) == dict: info = yi['info']   # information received which tells N is composite
                elif yi.degree % k: info = 'Sqrt step failed with deg(Fki) mod k != 0'  # save as above
                else: yi = yi.monic()
                Fki.append(yi)
                y = y*y % Fk
                if type(y) == dict: info = yi['info']   # information received which tells N is composite
        else: info = y['info']                          # information received which tells N is composite
        if info:               # if 'info' is not empty, N is composite
            if pr: print(info)
            if pr: print('COMPOSITE')
            return False
        y = 1
        for yi in Fki: y *= yi
        if y.coef != Fk.coef:
            if pr: print('Sqrt step failed with prod(Fki) != Fk')
            if pr: print('COMPOSITE')
            return False
        ## The additional code ends here
        fk //= Fk
        if k % 2 == 0: S += Fk.degree // k
    if fk.degree > 1:
        if pr: print('Factorization Step failed: fd(x) != 1')
        if pr: print('COMPOSITE')
        return False
    if (S%2, Jacobi(Delta,N)) not in [(0,1),(1,-1)]:
        if pr: print('Jacobi step failed.')
        if pr: print('COMPOSITE')
        return False
    if pr: print('Strong Frobenius probable prime test passed.')
    if pr: print('PRIME (probable)')
    return True

def Quadratic_Frobenius(N,b,c,pr=True):
    assert N > 1, 'Tested integer N should be NO LESS THAN 2'
    factor = 1
    ## Step 1: trial division
    for m in [2,3,5,7,11,13,17,19,23,29]:
        if m * m > N and factor == 1:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0: factor = m
    ## Step 2: 检验完全平方
    if powrt_int(N,2) > 0: factor = powrt_int(N,2)    # perfect square test
    ## 寻找其他可能的非平凡因数
    if 1 < gcd(b*b+4*c,N) < N: factor = gcd(b*b+4*c,N)
    if 1 < gcd(c,N) < N: factor = gcd(c,N)
    if factor > 1:
        if pr: print('Factor found: %d'%factor)
        if pr: print('COMPOSITE')
        return False
    assert Jacobi(b*b+4*c,N) == -Jacobi(-c,N) == -1, 'Condition for Jacobi symbol not satisfied.'
    x = Polynomial(modulo=N)
    f = x*x - b*x - c
    ## Step 3: 计算 x^{(N+1)/2} mod f(x)
    y = poly_pow_mod(x,N//2+1,f)
    if type(y) == dict: info = y['info']      # information received which tells N is composite
    elif y.degree > 0: info = 'Quadratic Frobenius test failed with x^((N+1)/2) not in Z/NZ.'    # Step 3 failed
    else: info = ''    # N is not found to be composite yet, so 'info' is empty
    if info:           # if 'info' is not empty, N is composite
        if pr: print(info)
        if pr: print('COMPOSITE')
        return False
    ## Step 4: 计算 x^(N+1)
    y = y.coef[0] if y.coef else 0
    if (y*y+c) % N > 0:            # Step 4 failed
        if pr: print('Quadratic Frobenius test failed with x^(N+1) != -c (mod N,f(x)).')
        if pr: print('COMPOSITE')
        return False
    ## Step 5,6: 
    r,s = 1,N*N-1
    while s % 2 == 0:
        r += 1
        s //= 2
    y = poly_pow_mod(x,s,f)
    if type(y) == Polynomial:
        if y.coef == {0:1}:       # x^s == 1 (mod N,f(x))
            if pr: print('Quadratic Frobenius test passed with x^s == 1 (mod N,f(x)).')
            if pr: print('PRIME (probable)')
            return True
        while not info:
            if y.coef == {0:N-1}:     # x^(s*2^i) == -1 (mod N,f(x))
                if pr: print('Quadratic Frobenius test passed with x^(s*2^i) == -1 (mod N,f(x)).')
                if pr: print('PRIME (probable)')
                return True
            r -= 1
            if r == 1: break
            y = y*y % f
            if type(y) == dict: info = y['info']       # information received which tells N is composite
        if r == 1: info = 'Quadratic Frobenius test failed with x^s != 1 and x^(s*2^i) != -1 (mod N,f(x)).'    # Step 5 failed
    else: info = y['info']          # information received which tells N is composite
    if pr: print(info)              # 'info' is not empty if the algorithm goes here, so N is composite
    if pr: print('COMPOSITE')
    return False

def Random_Quadratic_Frobenius(N,pr=True):
    factor = 1
    for m in [2,3,5,7,11,13,17,19]:            # trial division
        if m * m > N and factor == 1:
            if pr: print('Factor not found.')
            if pr: print('PRIME')
            return True
        if N % m == 0: factor = m
    if powrt_int(N,2): factor = powrt_int(N,2)    # perfect square test
    if factor > 1:
        if pr: print('Factor found: %d'%factor)
        if pr: print('COMPOSITE')
        return False
    while True:                                   # randomly choose b,c
        b,c = [random.randint(1,N-1) for _ in range(2)]
        if Jacobi(b*b+4*c,N) == -Jacobi(-c,N) == -1: break
    if pr: print('parameters b,c = %d,%d'%(b,c))
    return Quadratic_Frobenius(N,1,c,pr=pr)


##### 17. AKS #####

def AKS(N:int,pr=True):
    assert N > 1,'Tested integer should be NO LESS THAN 2'
    L = math.log2(N)
    ## Step 1
    b = 2
    while b <= L+1:
        a = powrt_int(N,b)
        if a > 0:
            if pr: print('Factor found: %d'%a)
            if pr: print('COMPOSITE')
            return False
        b += 1
    ## Step 2,3
    r = 2
    while r < N:
        if gcd(r,N) > 1:
            if pr: print('Factor found: %d'%r)
            if pr: print('COMPOSITE')
            return False
        if r > L*L:
            y,m = N%r,1
            while y != 1:
                y = y*N % r
                m += 1
            if m > L*L: break
        r += 1
    # Step 4
    if r == N:
        if pr: print('No factor found.')
        if pr: print('PRIME')
        return True
    ## Step 5
    x = Polynomial(modulo=N)
    p0 = Polynomial({0:-1,r:1},N)
    xN = poly_pow_mod(x,N,p0)
    phi = Euler_phi(r)
    max_a = math.floor(math.sqrt(phi)*L)
    if pr: print('r = %d,  phi(r) = %d'%(r,phi))
    if pr: print('max(a) = %d'%max_a)
    if pr: print('checking if for all a, (x+a)^N == x^N + a (mod N,x^r-1) ...')
    for a in range(1,max_a+1):
        p1 = poly_pow_mod(x+a,N,p0)
        if p1.coef != (xN+a).coef:
            if pr: print('AKS test failed with (x+a)^N != x^N + a (mod N,x^r-1), a,r = %d'%a)
            if pr: print('COMPOSITE')
            return False
        #print('a = %d, polynomial test passed'%a)
    ## Step 6
    if pr: print('AKS test passed.')
    if pr: print('PRIME')
    return True
