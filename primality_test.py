import math
import random

from .utils import (
    HINT_LARGER_THAN_1, HINT_NO_FACTOR, HINT_FACTOR,
    gcd, exgcd, power_mod, powrt_int, Jacobi,
    trial_div_factorize, Pollard_rho, FR_factorize, prob_prime_factorize,
    Cornacchia_4N, Lucas_sequence, V2_sequence, V3_sequence, Euler_phi,
    trial_div, check_sqrt
)
from .elliptic_curve import reduced_class_polynomial, complex_multiplication, Point, Elliptic_curve
from .cyclotomic import calculate_jpq, calculate_j2q, calculate_Jpq, calculate_J2q, find_apr_h
from .utils_apr import trial_division_30, apr_v, apr_e, generate_apr_f
from .polynomial import Polynomial, poly_pow_mod, gcd_poly, poly_nest


##### 0. Trial division #####

def Trial_division(N:int, pr:bool=True):
    assert N > 1, HINT_LARGER_THAN_1
    m = 2
    while m*m <= N:
        if N % m == 0:
            if pr: print(HINT_FACTOR, m)
            return False
        m += 1
    if pr: print(HINT_NO_FACTOR)
    return True

def Trial_division_30(N:int, pr:bool=True):
    assert N > 1, HINT_LARGER_THAN_1
    factor = 1
    for m in [2,3,5]:
        if m*m > N or factor > 1: break
        if N % m == 0: factor = m
    A,i,m = 0,1,7
    B = [1,7,11,13,17,19,23,29]
    while m*m <= N and factor == 1:
        if N % m == 0: factor = m
        i = (i+1) % 8
        if i == 0: A += 30
        m = A + B[i]
    if factor > 1:
        if pr: print(HINT_FACTOR, factor)
        return False
    if pr: print(HINT_NO_FACTOR)
    return True

def Trial_division_lim(N:int, pr:bool=True):
    assert N > 1, HINT_LARGER_THAN_1
    primes = []
    A,i,m = 2,0,2
    B = [1]
    k = r = 1
    while m*m <= N:
        m_prime = True
        for p in primes:
            if p*p > m: break
            if m % p == 0: m_prime = False
        if m_prime:
            primes.append(m)
            if N % m == 0:
                if pr: print(HINT_FACTOR, m)
                return False
        m = A*k + B[i]
        i += 1
        if i >= len(B):
            k += 1
            i = 0
        if len(primes) > r and k >= primes[r]:
            q = primes[r]
            new_B = []
            for k in range(q):
                for b in B:
                    if (A*k+b) % q > 0: new_B.append(A*k+b)
            A *= q
            r += 1
            B,k,i = new_B, 1, 0
    if pr: print(HINT_NO_FACTOR)
    return True


##### 1. Fermat #####

def fermat_core(N:int, a:int, pr:bool):
    if gcd(a,N) > 1:
        if pr: print(HINT_FACTOR, gcd(a,N))
        return False
    if power_mod(a, N-1, N) != 1:
        if pr: print('COMPOSITE : Fermat test failed with a =', a)
        return False
    return True

@ trial_div([2])
def Fermat(N:int, test_time:int=100, pr:bool=True):
    for _ in range(test_time):
        a = random.randint(2, N-2)
        if not fermat_core(N, a, pr): return False
    if pr: print(f'PRIME (probable) : Fermat test passed with {test_time} repeated tests.')
    return True


##### 2. Miller-Rabin #####

def miller_rabin_core(N:int, a:int, pr:bool):
    assert 1 < a < N-1
    if gcd(a,N) > 1:
        if pr: print(HINT_FACTOR, gcd(a,N))
        return False
    r,s = 0, N-1
    while s % 2 == 0:
        r += 1
        s //= 2
    y = power_mod(a, s, N)
    if y == 1: return True
    while True:
        if y == N-1: return True
        r -= 1
        if r == 0: break
        y = y*y % N
    if pr: print('COMPOSITE : Miller-Rabin test failed with a =', a)
    return False

@ trial_div([2])
def Miller_Rabin(N:int, test_time:int=100, pr:bool=True):
    for _ in range(test_time):
        a = random.randint(2, N-2)
        if not miller_rabin_core(N, a, pr): return False
    if pr: print(f'PRIME (probable) : Miller-Rabin test passed with {test_time} repeated tests.')
    return True


##### 3. Solovay-Strassen #####

def solovay_strassen_core(N:int, a:int, pr:bool):
    if gcd(a,N) > 1:
        if pr: print(HINT_FACTOR, gcd(a,N))
        return False
    if power_mod(a, N//2, N) != Jacobi(a,N) % N:
        if pr: print('COMPOSITE : Solovay-Strassen test failed with a =', a)
        return False
    return True

@ trial_div([2])
def Solovay_Strassen(N:int, test_time:int=100, pr:bool=True):
    for _ in range(test_time):
        a = random.randint(2, N-1)
        if not solovay_strassen_core(N, a, pr): return False
    if pr: print(f'PRIME (probable): Solovay-Strassen test passed with {test_time} repeated tests.')
    return True


##### 4. Lucas probable prime #####

def lucas_prob_prime_core(N:int, P:int, Q:int, strong:bool, extra:bool, pr:bool):
    assert strong or not extra, "Mode <strong> must be open when mode <extra> is open."
    D = P*P - 4*Q
    assert Q*D % N > 0, 'Tested integer N should NOT DEVIDE Q*D'
    if gcd(Q*D, N) > 1:
        if pr: print(HINT_FACTOR, gcd(Q*D, N))
        return False
    r,s = 0, N-Jacobi(D,N)
    while strong and s % 2 == 0:
        r += 1
        s //= 2
    U,V = Lucas_sequence(s,P,Q,N)
    if U == 0 and (not extra or V == 2 or V == N-2): return True
    while strong:
        if V == 0: return True
        r -= 1
        if r == 0: break
        U,V = U*V%N, D*U*U+V*V
        V = (V+V%2*N)//2 % N
    if extra: name = 'Extra strong Lucas probable prime test'
    elif strong: name = 'Strong Lucas probable prime test'
    else: name = 'Lucas probable prime test'
    if pr: print('COMPOSITE :', name, 'failed with (P,Q) =', (P,Q))
    return False

@ trial_div([2])
def Lucas_prob_prime(N:int, P:int, Q:int, pr:bool=True):
    if not lucas_prob_prime_core(N, P, Q, False, False, pr): return False
    if pr: print('PRIME (probable) : Lucas probable prime test passed with (P,Q) =', (P,Q))
    return True

@ trial_div([2])
def Lucas_prob_prime_random(N:int, test_time:int=100, pr:bool=True):
    test_time_left = test_time
    while test_time_left > 0:
        P, Q = random.randint(-200,200), random.randint(-200,200)
        try: ans = lucas_prob_prime_core(N, P, Q, False, False, pr)
        except AssertionError: continue
        if not ans: return False
        test_time_left -= 1
    if pr: print(f'PRIME (probable) : Lucas probable prime test passed with {test_time} repeated tests.')
    return True

@ trial_div([2])
def Strong_Lucas_prob_prime(N:int, P:int, Q:int, pr:bool=True):
    if not lucas_prob_prime_core(N, P, Q, True, False, pr): return False
    if pr: print('PRIME (probable) : Strong Lucas probable prime test passed with (P,Q) =', (P,Q))
    return True

@ trial_div([2])
def Strong_Lucas_prob_prime_random(N:int, test_time:int=100, pr:bool=True):
    test_time_left = test_time
    while test_time_left > 0:
        P, Q = random.randint(-200,200), random.randint(-200,200)
        try: ans = lucas_prob_prime_core(N, P, Q, True, False, pr)
        except AssertionError: continue
        if not ans: return False
        test_time_left -= 1
    if pr: print(f'PRIME (probable) : Strong Lucas probable prime test passed with {test_time} repeated tests.')
    return True

@ trial_div([2])
def Extra_strong_Lucas_prob_prime(N:int, b:int, pr:bool=True):
    if not lucas_prob_prime_core(N, b, 1, True, True, pr): return False
    if pr: print('PRIME (probable) : Extra strong Lucas probable prime test passed with (P,Q) =', (b,1))
    return True

@ trial_div([2])
def Extra_strong_Lucas_prob_prime_random(N:int, test_time:int=100, pr:bool=True):
    test_time_left = test_time
    while test_time_left > 0:
        b = random.randint(-200,200)
        try: ans = lucas_prob_prime_core(N, b, 1, True, True, pr)
        except AssertionError: continue
        if not ans: return False
        test_time_left -= 1
    if pr: print(f'PRIME (probable) : Extra strong Lucas probable prime test passed with {test_time} repeated tests.')
    return True


##### 5. Baillie-PSW #####

@ trial_div([2,3,5,7,11,13,17,19])    # trial division
@ check_sqrt                          # perfert square test
def Baillie_PSW(N:int, pr:bool=True):
    ## Miller-Rabin based on 2
    if not miller_rabin_core(N, 2, pr): return False
    if pr: print('Miller-Rabin test (basis 2) passed.')
    ## Find D
    D = 5
    P,Q = 1, (1-D)//4
    while Jacobi(D,N) > -1 or Q*D % N == 0:
        D = -(D+2) if D > 0 else -(D-2)
        P,Q = 1, (1-D)//4
    ## strong Lucas probable prime test
    if not lucas_prob_prime_core(N, P, Q, True, False, pr): return False
    V = Lucas_sequence(N+1, P, Q, N)[1]
    if pr:
        if V == 2*Q%N: print('PRIME (probable) : Strong Lucas probable prime test passed with (P,Q) =', (P,Q))
        else: print(f'COMPOSITE : Strong Lucas probable prime test failed with (P,Q) = {(P,Q)} ,  V_(N+1) != 2*Q (mod N)')
    return V == 2*Q%N


##### 6. Lucas #####

def lucas_core(N:int, Qs:set, a:int, pr:bool):
    for q in Qs:
        if power_mod(a, (N-1)//q, N) == 1:
            return False
    if pr: print('PRIME : Lucas test passed with a =', a)
    return True

@ trial_div([2])
def Lucas(N:int, test_time:int=100, pr:bool=True):
    Qs = trial_div_factorize(N-1)
    if pr: print(f'Prime factors of N-1: {Qs}')
    for _ in range(test_time):
        a = random.randint(2, N-2)
        if not fermat_core(N, a, pr): return False
        if lucas_core(N, Qs, a, pr): return True
    if pr: print(f'COMPOSITE (probable) : Lucas test failed with {test_time} repeated tests.')
    return False


##### 7. Pocklington #####

def pocklington_core(N:int, Qs:set, test_time:int, pr:bool):
    for q in Qs:
        found_a = False
        for _ in range(test_time):
            a = random.randint(2, N-2)
            if not fermat_core(N, a, pr): return False
            x = power_mod(a, (N-1)//q, N) - 1
            if gcd(x,N) == 1:
                found_a = True
                break
        if not found_a:
            if pr: print('COMPOSITE (probable) : Pocklington test failed with no basis found for q =', q)
            return False
    if pr: print('PRIME : Pocklington test passed.')
    return True

@ trial_div([2,3,5,7,11,13,17,19])
def Pocklington(N:int, known_factors:set=set(), test_time:int=100, pr:bool=True):
    R = N-1
    Qs = known_factors.copy()
    for q in Qs:
        assert Trial_division_30(q, False), f'Known factors should be prime, but {q} is not.'
        if R % q > 0:
            Qs.remove(q)
            print(f'Known factors should devide N-1, but {q} does not. It has been removed from known factors.')
        while R % q == 0: R //= q
    if (N-1)//R <= R:
        r = FR_factorize(R, balance=True)[0]
        R //= r
        Qs |= trial_div_factorize(r)
    if pr: print(f'Factorization of N-1: F = {(N-1)//R}, R = {R}')
    if pr: print('Prime factors of F:', Qs)
    return pocklington_core(N, Qs, test_time, pr)


##### 8. Proth (for P=A*2^n+1) & Pépin (for F_n=2^{2^n}+1) #####

def Proth(A:int, n:int, test_time:int=100, pr:bool=True):
    assert n > 0, 'Integer n should be POSITIVE.'
    assert A > 0 and A % 2 > 0 and A < 2**n, 'Integer A should be POSITIVE ODD and LESS than 2^n.'
    P = A * 2**n + 1
    if pr: print('Tested integer: P =', P)
    for _ in range(test_time):
        a = random.randint(2, P-1)
        if gcd(a,P) > 1:
            if pr: print(HINT_FACTOR, gcd(a,P))
            return False
        y = power_mod(a,A,P)
        for _ in range(n-1): y = y*y % P
        if y == P-1:
            if pr: print('PRIME : Proth test passed with a =', a)
            return True
        if y == 1 == -Jacobi(a,P) or y != 1:
            if pr: print('COMPOSITE : Fermat or Solovay-Strassen test failed with a =', a)
            return False
    if pr: print(f'COMPOSITE (probable) : Proth test failed with {test_time} repeated tests.')
    return False

def Pepin(n:int, pr:bool=True):
    assert n >= 0, 'Integer n should be NOT NEGATIVE.'
    m = 2**n
    Fn = 2**m + 1
    if pr: print('Tested integer: Fn =', Fn)
    if n == 0:
        if pr: print('PRIME')
        return True
    y = 3
    for _ in range(m-1): y = y*y % Fn
    if pr:
        if y == Fn-1: print('PRIME : 3^((Fn-1)/2) == -1 (mod Fn)')
        else: print('COMPOSITE : 3^((Fn-1)/2) != -1 (mod Fn)')
    return y == Fn-1


##### 9. Factorizing Lucas probable prime #####

@ trial_div([2,3,5,7,11,13,17,19])
@ check_sqrt
def Factorizing_Lucas_prob_prime_v0(N:int, test_time:int=100, pr:bool=True):
    Qs = trial_div_factorize(N+1)
    if pr: print('Prime factors of N+1:', Qs)
    D = 5
    while test_time:
        P,Q = 1,(1-D)//4
        if Jacobi(D,N) > -1 or Q*D % N == 0:
            D = -(D+2) if D > 0 else -(D-2)
            continue
        if gcd(Q*D, N) > 1:
            if pr: print(HINT_FACTOR, gcd(Q*D, N))
            return False
        U,V = Lucas_sequence(N+1, P, Q, N)
        if U > 0 or V != 2*Q%N:
            if pr: print('COMPOSITE : Lucas probable prime test failed with (P,Q) =', (P,Q))
            return False
        q_test = True
        for q in Qs:
            if Lucas_sequence((N+1)//q, P, Q, N)[0] == 0:
                q_test = False
                break
        if q_test:
            if pr: print(f'PRIME : Factorizing Lucas probable prime test (v0) passed with (P,Q) =', (P,Q))
            return True
        D = -(D+2) if D > 0 else -(D-2)
        test_time -= 1
    if pr: print(f'COMPOSITE (probable) : Factorizing Lucas probable prime test (v0) failed with {test_time} repeated tests.')
    return False

@ trial_div([2,3,5,7,11,13,17,19])
@ check_sqrt
def Factorizing_Lucas_prob_prime_v1(N:int, known_factors:set=set(), test_time:int=100, pr:bool=True):
    R = N+1
    Qs = known_factors.copy()
    for q in Qs:
        assert Trial_division_30(q, False), f'Known factors should be prime, but {q} is not.'
        if R % q > 0:
            Qs.remove(q)
            print(f'Known factors should devide N+1, but {q} does not. It has been removed from known factors.')
        while R % q == 0: R //= q
    F = (N+1) // R
    if (F-1)**2 <= N:
        for _ in range(10):
            q = FR_factorize(R, balance=True)[0]
            if (F*q-1)**2 > N: break
        if (F*q-1)**2 <= N: q = R
        F,R = F*q, R//q
        Qs |= trial_div_factorize(q)
    if pr: print(f'Factorization of N-1: F = {F}, R = {R}')
    if pr: print('Prime factors of F:', Qs)
    for q in Qs:
        found_PQ = False
        D = 5
        while test_time:
            P,Q = 1,(1-D)//4
            if Jacobi(D,N) > -1 or Q*D % N == 0:
                D = -(D+2) if D > 0 else -(D-2)
                continue
            if gcd(Q*D, N) > 1:
                if pr: print(HINT_FACTOR, gcd(Q*D, N))
                print(False)
            U,V = Lucas_sequence(N+1, P, Q, N)
            if U > 0 or V != 2*Q%N:
                if pr: print('COMPOSITE : Lucas probable prime test failed with (P,Q) =', (P,Q))
                return False
            u = Lucas_sequence((N+1)//q, P, Q, N)[0]
            if gcd(u,N) == 1:
                found_PQ = True
                break
            D = -(D+2) if D > 0 else -(D-2)
            test_time -= 1
        if not found_PQ:
            if pr: print('COMPOSITE (probable) : Factorizing Lucas probable prime test failed with no (P,Q) found for q =', q)
            return False
    if pr: print('PRIME : Factorizing Lucas probable prime test passed.')
    return True


##### 10. Lucas-Lehmer (for M_p=2^p-1) & Lucas-Lehmer-Riesel (for N=A*2^n-1) #####

def Lucas_Lehmer(p:int, pr:bool=True):
    assert p > 1, 'Integer p should be larger than 1.'
    Mp = 2**p - 1
    if pr: print('Tested integer: Mp =', Mp)
    if p == 2:
        if pr: print('PRIME')
        return True
    # trial test for p
    factor = 1
    for m in [2,3,5]:
        if m*m > p or factor > 1: break
        if p % m == 0: factor = m
    A,i,m = 0,1,7
    B = [1,7,11,13,17,19,23,29]
    while m*m <= p and factor == 1:
        if p % m == 0: factor = m
        i = (i+1) % 8
        if i == 0: A += 30
        m = A + B[i]
    if factor > 1:
        if pr: print(HINT_FACTOR, 2**factor-1)
        return False
    # Lucas_Lehmer for Mp
    S = 4 % Mp
    for _ in range(p-2): S = (S*S-2) % Mp
    if pr:
        if S > 0: print('COMPOSITE : Sp != 0 (mod Mp)')
        else: print('PRIME : Sp == 0 (mod Mp)')
    return S == 0

def Lucas_Lehmer_Riesel(A:int, n:int, pr:bool=True):
    assert A % 2 > 0 and 2 < A < 2**n, 'Integer A should be ODD and within [3, 2**n-1].'
    assert n > 1, 'Integer n should be larger than 1.'
    N = A * 2**n - 1
    if pr: print('Tested integer: N =', N)
    if A % 3 > 0: P = 4
    else:
        P = 5
        while Jacobi(P-2,N) != 1 or Jacobi(P+2,N) != -1:
            if P == 5: P = 8
            elif P == 9: P = 11
            else: P += 1
    S = Lucas_sequence(A,P,1,N)[1]
    for _ in range(n-2): S = (S*S-2) % N
    if pr:
        if S > 0: print('COMPOSITE : Sn != 0 (mod N)')
        else: print('PRIME : Sn == 0 (mod N)')
    return S == 0


##### 11. ECPP: Goldwasser-Kilian #####

@ trial_div([2])
def Goldwasser_Kilian(N:int, test_time=50, pr:bool=True):
    assert N > 1, HINT_LARGER_THAN_1
    test_number = 1
    while test_number <= test_time:
        if pr: print(f'Test No.{test_number}')
        ## Step 1: randomly choose a,b,x,y
        a,x,y = [random.randint(0,N-1) for _ in range(3)]
        b = (y*y - x**3 - a*x) % N
        if (4*a**3 + 27*b*b) % N == 0: continue
        P = Point(x,y)
        if pr: print(f'a,b = {a},{b}   P = {P}')
        test_number += 1
        ## Step 2: Schoof's algorithm to get order m
        ep = Elliptic_curve(N,a,b)
        if pr: print('running Schoof\'s algorithm ...')
        try: m = ep.schoof()
        except ArithmeticError as info:
            if pr: print('\nCOMPOSITE :', info)
            return False
        except ValueError as x:
            if pr: print('\n'+HINT_FACTOR, x)
            return False
        if pr: print('Order of group: m =', m)
        ## Step 3: Find factor s of m, and factorize s
        for _ in range(10):
            s = FR_factorize(m,balance=True)[0]
            if s > (N**(1/4)+1)**2: break
        if s <= (N**(1/4)+1)**2: s = m
        if pr: print('factor of m: s =', s)
        Q1,Q2 = prob_prime_factorize(s)
        if pr: print('factors of s:', Q1, Q2)
        ## Steo 4,5,6: Compute mP and (m/q)P
        try: Pm = ep.times(m,P)
        except ArithmeticError as info:
            if pr: print('\nCOMPOSITE :', info)
            return False
        except ValueError as x:
            if pr: print('\n'+HINT_FACTOR, x)
            return False
        if Pm != 0:
            if pr: print('\nCOMPOSITE : Group order test failed: mP != 0')
            return False
        if pr: print('mP = 0')
        q_test = True
        for q in Q1|Q2:
            try: Pm = ep.times(m//q, P)
            except ArithmeticError as info:
                if pr: print('\nCOMPOSITE :', info)
                return False
            except ValueError as x:
                if pr: print('\n'+HINT_FACTOR, x)
                return False
            if pr: print(f'(m/{q})P = {Pm}')
            if Pm == 0:
                if pr: print('Let\'s try again ...\n')
                q_test = False
                break
        if not q_test: continue
        ## Step 7: Recursion primality test for all q in Q2
        q_test = True
        for q in Q2:
            if pr: print(f'\nRecursion: ECPP for q = {q} ...')
            if not Goldwasser_Kilian(q, pr=pr):
                if pr: print(f'\nq = {q} is composite, try again ...\n')
                q_test = False
                break
            elif pr: print(f'\nq = {q} is prime.')
        if not q_test: continue
        if pr: print(f'\nPRIME : ECPP success with a,b,x,y = {a},{b},{x},{y}')
        return True
    ## Step 8: Sufficient test to tell N is composite
    if pr: print('\nCOMPOSITE (probable) : Elliptic curve with proper order not found.')
    return False


##### 12. ECPP: Atkin-Morain #####

@ trial_div([2,3,5,7])
def Atkin_Morain(N:int, max_D=500, pr:bool=True):
    assert N > 1, HINT_LARGER_THAN_1
    square_factors = [i*i for i in range(2,int(math.sqrt(max_D))+1)]
    test_number = 1
    for D in random.sample(range(1,max_D+1),max_D):
        ## Step 1: choose D, and solve equation 4N == u^2 + Dv^2
        D = random.randint(1, max_D)
        if 1 < gcd(D, N) < N:
            if pr: print(HINT_FACTOR, gcd(D, N))
            return False
        if any(D % s == 0 for s in square_factors): continue
        if Jacobi(-D,N) < 1: continue
        if pr: print(f'Test No.{test_number}:  D = {D}')
        test_number += 1
        try: RCP = reduced_class_polynomial(D)
        except Exception as info:
            if pr: print(info,'\n')
            continue
        RCP.modulo = N
        try: uvs = Cornacchia_4N(N,D)
        except ArithmeticError as info:
            if pr: print('\nCOMPOSITE :', info)
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
        if pr: print(f'u = {u}, probable values of order m: {ms}')
        ## Step 2: complex multiplication method to get a,b, and randomly choose a point P on the curve
        a,b = complex_multiplication(N, D, RCP, u%2)
        if pr: print(f'a = {a}, b = {b}')
        ep = Elliptic_curve(N,a,b)
        try: P = ep.random_point()
        except ArithmeticError as info:
            if pr: print('\nCOMPOSITE :', info)
            return False
        except ValueError as x:
            if pr: print('\n'+HINT_FACTOR, x)
            return False
        if pr: print(f'randomly choose a point: P(x,y) = {P}')
        ## Step 3: find m such that mP == 0
        for m in ms+[0]:
            try: Pm = ep.times(m,P)
            except ArithmeticError as info:
                if pr: print('\nCOMPOSITE :', info)
                return False
            except ValueError as x:
                if pr: print('\n'+HINT_FACTOR, x)
                return False
            if Pm == 0: break                        # mP == 0
        if m == 0:                                               # all mP != 0, return composite
            if pr: print('\nCOMPOSITE : Group order test failed: mP != 0')
            return False
        if pr: print('Order of group: m =', m)
        ## Step 4: Find factor s of m, and factorize s
        for _ in range(10):
            s = FR_factorize(m,balance=True)[0]
            if s > (N**(1/4)+1)**2: break
        if s <= (N**(1/4)+1)**2: s = m
        if pr: print(f'factor of m: s =', s)
        Q1,Q2 = prob_prime_factorize(s)
        if pr: print('factors of s:', Q1, Q2)
        ## Steo 5,6: Compute (m/q)P
        q_test = True
        for q in Q1|Q2:
            try: Pm = ep.times(m//q, P)
            except ArithmeticError as info:
                if pr: print('\nCOMPOSITE :', info)
                return False
            except ValueError as x:
                if pr: print('\n'+HINT_FACTOR, x)
                return False
            if pr: print(f'(m/{q})P = {Pm}')
            if Pm == 0:                          # (m/q)P == 0, back to step 1 and choose another D
                if pr: print('Let\'s try again ...\n')
                q_test = False
                break
        if not q_test: continue
        ## Step 7: Recursion primality test for all q in Q2
        q_test = True
        for q in Q2:
            if pr: print(f'\nRecursion: ECPP for q = {q} ...')
            if not Atkin_Morain(q, pr=pr):          # one of the probable prime factors is said to be composite
                if pr: print(f'\nq = {q} is composite.')
                q_test = False
                break
            elif pr: print(f'\nq = {q} is prime, try again ...\n')
        if not q_test: continue
        if pr: print(f'\nPRIME : ECPP success with m,a,b,x,y = {m},{a},{b},{P.x},{P.y}')    # All probable prime factors are definitely prime, so N is prime
        return True
    ## Step 8: Sufficient test to tell N is composite
    if pr: print('\nCOMPOSITE (probable) : Elliptic curve with proper order not found.')
    return False


##### 13. ECPP: Tsumura #####

def Tsumura(k:int, n:int, pr:bool=True):
    assert k >= 2, 'Integer k should be larger than 1.'
    assert n >= 1 and n % 2 > 0, 'Integer n should be POSITIVE and ODD.'
    N = 2**k * n - 1
    M = (N**0.25+1)**2
    assert n < N/M or n > M, 'Integer n not in applicable range.'
    if pr: print('Tested integer: N =', N)
    while True:
        x = random.randint(1, N-1)
        jacobi = Jacobi(x,N)
        if jacobi == -1: break
        if jacobi == 0:
            if pr: print(HINT_FACTOR, gcd(x,N))
            return False
    if pr: print('randomly choose x =', x)
    while True:
        while True:
            y = random.randint(1, N-1)
            m = (x**3-y*y) % N
            if m == 0: continue
            jacobi = Jacobi(m,N)
            if jacobi == 1: break
            if jacobi == 0:
                if pr: print(HINT_FACTOR, gcd(m,N))
                return False
        m = m * exgcd(x,N)[1] % N
        ep = Elliptic_curve(N,-m,0)
        P = Point(x,y)
        if pr: print(f'randomly choose y = {y}\nm = {m}')
        if n < N/M:
            if pr: print('Range of n:  n < N/(N^(1/4)+1)^2')
            try: R = ep.times(n,P)
            except ArithmeticError as info:
                if pr: print('COMPOSITE :', info)
                return False
            except ValueError as x:
                if pr: print(HINT_FACTOR, x)
                return False
            X = R.x
            for i in range(1, k):
                S = X**3 - m*X
                if S % N == 0:
                    if pr: print(f'COMPOSITE : Order test failed: (2^{i})nP == 0')
                    return False
                if gcd(S,N) > 1:
                    if pr: print(HINT_FACTOR, gcd(S,N))
                    return False
                X = (X*X+m)**2 * exgcd(4*S,N)[1] % N
            S = (X**3 - m*X) % N
            if pr:
                if S > 0: print(f'COMPOSITE : Order test failed: (2^{k})nP != 0')
                else: print(f'PRIME : ECPP success with m,x,y = {m},{x},{y}')
            return S == 0
        if pr: print('Range of n:  n > (N^(1/4)+1)^2')
        R = P
        for _ in range(k):
            try: R = ep.plus(R,R)
            except ArithmeticError as info:
                if pr: print('COMPOSITE :', info)
                return False
            except ValueError as x:
                if pr: print(HINT_FACTOR, x)
                return False
        if pr: print('(2^k)P =', R)
        if R == 0:
            if pr: print('try another y ...')
            continue
        Qs = trial_div_factorize(n)
        if pr: print('prime factors of n:', Qs)
        q_test = True
        for q in Qs:
            try: Rn = ep.times(n//q, R)
            except ArithmeticError as info:
                if pr: print('COMPOSITE :', info)
                return False
            except ValueError as x:
                if pr: print(HINT_FACTOR, x)
                return False
            if pr: print(f'(2^k)(n/{q})P = {Rn}')
            if Rn == 0:
                if pr: print('try another y ...')
                q_test = False
                break
        if not q_test: continue
        try: Rn = ep.times(n,R)
        except ArithmeticError as info:
            if pr: print('COMPOSITE :', info)
            return False
        except ValueError as x:
            if pr: print(HINT_FACTOR, x)
            return False
        if pr:
            if Rn == 0: print(f'PRIME : ECPP success with m,x,y = {m},{x},{y}')
            else: print('COMPOSITE : Order test failed: (2^k)nP != 0')
        return Rn == 0


##### 14. High-order field & cyclotomic polynomial #####

@ trial_div([2,3,5,7,11,13,17,19])
def Quartic_field(N:int, max_D=120, max_C=120, test_time:int=100, pr:bool=True):
    if pr: print('Phi = N*N + 1 =', N*N+1)
    F4, R4 = 1, N*N//2+1
    while (F4+1)**3 <= N*N:
        B4 = (N**(2/3)+1) / F4
        s = Pollard_rho(R4, test_time=10, loop_time=max(int(math.sqrt(B4)),2**16))
        if s == 1: break
        R4 //= s
        r = gcd(s,R4)
        while r > 1:
            s *= r
            R4 //= r
            r = gcd(s,R4)
        F4 *= s
    if pr: print(f'F4 = {F4}, R4 = {R4}')
    Q4 = trial_div_factorize(F4)
    if R4 > 1: Q4.add(R4)
    if pr: print('All factors to be tested:', Q4)
    D = 5
    while abs(D) < max_D:
        if 1 < gcd(D,N) < N:
            if pr: print(HINT_FACTOR, gcd(D,N))
            return False
        if D % N == 0 or Jacobi(D,N) == 1:
            D = -(D+2) if D > 0 else -(D-2)
            continue
        for C in range(max_C):
            if 1 < gcd(C*C-16*D, N) < N:
                if pr: print(HINT_FACTOR, gcd(C*C-16*D, N))
                return False
            if (C*C-16*D) % N == 0 or Jacobi(C*C-16*D, N) == 1: continue
            if pr: print(f'\nSet parameter: D,C = {D},{C}')
            q_test = True
            for q in Q4:
                param_found = False
                for _ in range(test_time):
                    H,K = [random.randint(0,N-1) for _ in [0,1]]
                    if (H*H-K*K*D) % N == 0: continue
                    if gcd(H*H-K*K*D, N) > 1:
                        if pr: print(HINT_FACTOR, gcd(H*H-K*K*D, N))
                        return False
                    p = (2*H*H + H*K*C + 2*K*K*D) % N
                    P1,P2,Q = 4*p % N, 4*(p*p-D) % N, (p*p-H*H*C-8*H*K*D-K*K*C*D+D) % N
                    if V2_sequence(N*N+1, P1, P2, Q, N)[1] % N > 0:
                        if pr: print(f'COMPOSITE : Quartic field test failed: V_(N^2+1) != 0 (mod N), with D,C,H,K = {D},{C},{H},{K}')
                        return False
                    V = V2_sequence((N*N+1)//q, P1, P2, Q, N)[1]
                    if gcd(V,N) == 1:
                        param_found = True
                        break
                if not param_found:
                    q_test = False
                    break
            if q_test:
                if pr: print(f'PRIME : Quartic field test passed with D,C = {D},{C}')
                return True
        D = -(D+2) if D > 0 else -(D-2)
    if pr: print('COMPOSITE (probable) : Proper parameters not found.')
    return False

@ trial_div([2,3,5,7,11,13,17,19])
def cubic_sextic_field(N:int, theta:int, max_G, test_time, pr):
    assert abs(theta) == 1, 'Parameter theta must be 1 or -1.'
    Phi = N*N + theta*N + 1
    if pr: print('Phi = N^2', '+' if theta>0 else '-', 'N + 1 =', Phi)
    F,R = 1,Phi
    while F**3 <= N*N:
        s = Pollard_rho(R, test_time=10, loop_time=max(int(R**(1/4)),1<<16))
        if s == 1: s = R
        R //= s
        r = gcd(s,R)
        while r > 1:
            s *= r
            R //= r
            r = gcd(s,R)
        F *= s
    if pr: print(f'F = {F}, R = {R}')
    Qs = trial_div_factorize(F)
    if pr: print('All factors to be tested:', Qs)
    for P,S,T in [(7,1,1),(19,7,1),(31,4,2),(61,1,3),(73,7,3),(97,19,1),(103,13,3),(151,19,3),(163,25,1),
                  (181,7,5),(211,13,5),(223,28,2),(229,22,4),(307,16,6),(331,1,7),(349,37,1),(373,13,7),
                  (397,34,4),(409,31,5),(421,19,7),(439,28,6),(457,10,8),(487,25,7)]:
        if N == P:
            if pr: print('PRIME : Contained in prime list.')
            return True
        if gcd(P,N) > 1:
            if pr: print(HINT_FACTOR, gcd(P,N))
            return False
        if power_mod(N%P, (P-1)//3, P) == 1: continue
        for G in range(1,max_G):
            if gcd(G,N) > 1:
                if pr: print(HINT_FACTOR, gcd(G,N))
                return False
            if G % N == 0 or Jacobi(G,N) == -theta: continue
            if pr: print(f'\nSet parameter: P,S,G = {P},{S},{G}')
            q_test = True
            for q in Qs:
                param_found = False
                for _ in range(test_time):
                    H,K,L = [random.randint(0,N-1) for _ in range(3)]
                    a,b = 3*P, P*S
                    m1,m2,m3 = (H*H+2*b*L*K-G) % N, (b*L*L+2*H*K+2*a*K*L) % N, (K*K+a*L*L+2*H*L) % N
                    d1,d2,d3 = ((m1+a*m3)**2 - m2*(a*m2+b*m3)) % N, (b*m3*m3-m1*m2) % N, (m2*m2-m1*m3-a*m3*m3) % N
                    r = (d1*m1 + b*d2*m3 + b*d3*m2) % N
                    A,B,C = (4*G*d1+2*r) % N, 4*G*d2 % N, 4*G*d3 % N
                    P1,P2,P3,Q = (3*A+2*a*C) % N, ((A+a*C)*(3*A+a*C)-a*B*B-3*b*B*C) % N, (A*(A+a*C)**2+b*B**3+b*b*C**3-3*b*A*B*C-a*A*B*B-a*b*B*C*C) % N, r*r % N
                    Delta = (P1*P2*(P1*P2+18*P3) - 4*P2**3 - P3*(4*P1**3+27*P3))%N
                    E = ((P3+4*Q*P1)**2 - 4*Q*(P2+4*Q)**2) % N
                    if (Delta*E*Q) % N == 0: continue
                    if gcd(Delta*E*Q, N) > 1:
                        if pr: print(HINT_FACTOR, gcd(Delta*E*Q, N))
                        return False
                    _,v1,v2 = V3_sequence(Phi,P1,P2,P3,Q,N)
                    if gcd(v1,v2) % N > 0:
                        if pr: print('COMPOSITE :', 'Cubic' if theta>0 else 'Sextic', f'field test failed with P,S,G,H,K,L = {P},{S},{G},{H},{K},{L}')
                        return False
                    _,v1,v2 = V3_sequence(Phi//q, P1, P2, P3, Q, N)
                    if gcd(v1,v2,N) == 1:
                        param_found = True
                        break
                if not param_found:
                    q_test = False
                    break
            if q_test:
                if pr: print('PRIME :', 'Cubic' if theta>0 else 'Sextic', f'field test success with P,S,G = {P},{S},{G}')
                return True
    if pr: print('COMPOSITE (probable) : Parameters not found.')
    return False

def Cubic_field(N:int, max_G=120, test_time:int=100, pr:bool=True):
    return cubic_sextic_field(N,theta=1, max_G=max_G, test_time=test_time, pr=pr)

def Sextic_field(N:int, max_G=120, test_time:int=100, pr:bool=True):
    return cubic_sextic_field(N, theta=-1, max_G=max_G, test_time=test_time, pr=pr)


##### 15. APR-CL test #####

@ trial_div([2,3,5,7,11,13,17,19])
def APR_CL(N:int, pr:bool=True):
    ## Step 1
    for t in [2,6,28,20,12,24,40,48,36,60,72,120,300,180,360,420,720,840,1080,1200,
              1260,1680,2520,7920,5040,12600,15120,25200,55440,113400,110880,221760,504000,720720,
              1441440,4324320,24504480,73513440,367567200,1396755360,6983776800]:      # A list of parameter t, which supports tests for at most 6011 digits
        s = apr_e(t)
        if s * s > N: break
    if pr: print(f'Step 1\nt,s = {t},{s}')
    if gcd(t*s, N) > 1:
        if pr: print(HINT_FACTOR, gcd(t*s,N))
        return False
    Euclidean = trial_div_factorize(s)                     # list all Euclidean primes q
    if pr: print('Euclidean primes:', Euclidean)
    ## Step 2
    if pr: print('\nStep 2')
    js, js_ast, js_tag = {}, {}, {}
    for q in Euclidean:
        if q == 2: continue
        f = generate_apr_f(q)                     # Step 2.1
        for p in trial_div_factorize(q-1):        # Step 2.2
            k = apr_v(p,q-1)                        # Step 2.2.1
            if p**k > 2: js[(p,q)] = calculate_jpq(p,k,q,f)    # Step 2.2.2
            if p == 2 and k > 2: js_ast[q],js_tag[q] = calculate_j2q(k,q,js[(2,q)],f)   # Step 2.2.3
        if pr: print(f'j(p,q) for q = {q} are calculated.')
    ## Stage B: loop of Step 3,4,5
    if pr: print('\nStep 3,4,5')
    for p in trial_div_factorize(t):
        lbd = p > 2 and power_mod(N, p-1, p*p) != 1       # Step 3
        for q in Euclidean:                               # Step 4
            if (q-1) % p > 0: continue
            k = apr_v(p,q-1)                              # Step 4.1
            uk,vk = N//p**k, N%p**k
            if p > 2: J0,Jv = calculate_Jpq(p, k, vk, js[(p,q)], N)                 # Step 4.2
            elif k == 1: J0,Jv = q,1
            elif k == 2: J0,Jv = js[(p,q)]**2*q, 1 if vk==1 else js[(p,q)]**2
            else: J0,Jv = calculate_J2q(k, vk, js_ast[q], js_tag[q],N)
            z = power_mod(J0,uk,N) * Jv                  # Step 4.3
            h = find_apr_h(p,k,z,N)                       # find h
            if h == p**k:
                if pr: print(f'\nCOMPOSITE : APR-CL test failed at Step 4.3 with p,q = {p},{q}')
                return False
            if h % p > 0 and (p > 2 or k == 1 and N % 4 == 1): lbd = True            # Step 4.4
            if h % 2 > 0 and p == 2 and k > 1 and not lbd:                           # Step 4.5
                if power_mod(q,N//2,N) < N-1:
                    if pr: print(f'\nCOMPOSITE : APR-CL test failed at Step 4.5 with p,q = {p},{q}')
                    return False
                lbd = True
        if lbd and pr: print(f'lambda_(p={p}) is True')
        if not lbd:                                           # Step 5
            if pr: print(f'lambda_(p={p}) is False ...')
            q = 2*p+1 if p>2 else 5                        # Step 5.1
            try_time = 0
            while try_time < len(Euclidean) + 100:              # max attempt times
                if trial_division_30(q):                     # simply trial division
                    if q not in Euclidean and power_mod(N, (q-1)//p, q) != 1: break
                    try_time += 1
                q += 2*p if p>2 else (4 if q%6==1 else (8 if N%4==3 else 2))
            if try_time == len(Euclidean) + 100:               # q' not found within max attempt times, algorithm fails
                if pr: print('\nUnable to test primality as q\' not found.')
                return
            if pr: print(f'prime q\' = {q} is found')
            if N % q == 0:                  # Step 5.2
                if pr: print(HINT_FACTOR, q)
                return False
            k = 2 if p==2 and N%4==3 else 1       # Step 5.3
            uk,vk = N//p**k, N%p**k
            f = generate_apr_f(q)
            if p**k > 2: j = calculate_jpq(p,k,q,f)
            if p > 2: J0,Jv = calculate_Jpq(p,1,vk,j,N)
            elif k == 1: J0,Jv = q,1
            else: J0,Jv = j*j*q, j*j
            if pr: print(f'j(p,q\') for q\' = {q} are calculated.')
            z = power_mod(J0,uk,N) * Jv      # Step 5.4
            h = find_apr_h(p,k,z,N)
            if h % p == 0:
                if pr: print(f'\nCOMPOSITE : APR-CL test failed at Step 5.4 with p,q\' = {p},{q}')
                return False
            if k == 2 and power_mod(q,N//2,N) < N-1:       # Step 5.5
                if pr: print(f'\nCOMPOSITE : APR-CL test failed at Step 5.5 with p,q\' = {p},{q}')
                return False
            if pr: print(f'lambda_(p={p}) -> True')
    ## Stage C: Step 6,7
    if pr: print('\nStep 6,7')
    r = 1
    for _ in range(1,t+1):
        r = r*N % s
        if r == 1:
            if pr: print('\nPRIME : APR-CL test success.')
            return True
        if N % r == 0:
            if pr: print(HINT_FACTOR, r)
            return False


##### 16. Frobenius probable prime test #####

def frobenius_prob_prime_core(N:int, f:Polynomial, strong:bool, pr:bool=True):
    f.modulo = N
    d = f.degree
    Delta = f.discriminant()
    assert f.coef[0]*Delta % N > 0, 'Tested interger N should not devide f(0)*Delta'
    if gcd(f.coef[0]*Delta, N) > 1:
        if pr: print(HINT_FACTOR, gcd(f.coef[0]*Delta, N))
        return False
    x = Polynomial(modulo=N)
    fk = f.copy()
    Nk = 1
    S = 0
    for k in range(1,d+1):
        Nk *= N
        try:
            Fk = poly_pow_mod(x,Nk,fk) - x
            Fk = gcd_poly(Fk,fk)
            Fk = Fk.monic()
            y = poly_pow_mod(x,N,Fk)
            y = poly_nest(Fk,y) % Fk
        except ValueError as x:                               # information received which tells N is composite
            if pr: print(HINT_FACTOR, x)
            return False
        if y.coef:
            if pr: print('COMPOSITE : Frobenius step failed at k =', k)
            return False
        ## additional square-root-step for strong test
        if strong:
            r,s = 1,Nk-1
            while s % 2 == 0:
                r += 1
                s //= 2
            try:
                y = poly_pow_mod(x,s,Fk)
                yi = gcd_poly(y-1,Fk)               # compute Fk0(x)
                yi = yi.monic()
                Fki = [yi]
                if yi.degree % k > 0: raise AttributeError('Square-root step failed with deg(Fki) mod k != 0')
                while len(Fki) <= r:
                    yi = gcd_poly(y+1,Fk)
                    if yi.degree % k > 0: raise AttributeError('Square-root step failed with deg(Fki) mod k != 0')
                    yi = yi.monic()
                    Fki.append(yi)
                    y = y*y % Fk
            except AttributeError as info:               # if 'info' is not empty, N is composite
                if pr: print('COMPOSITE :', info)
                return False
            except ValueError as x:
                if pr: print(HINT_FACTOR, x)
                return False
            y = 1
            for yi in Fki: y *= yi
            if y.coef != Fk.coef:
                if pr: print('COMPOSITE : Square-root step failed with prod(Fki) != Fk')
                return False
        fk //= Fk
        if k % 2 == 0: S += Fk.degree // k
    if fk.degree > 1:
        if pr: print('COMPOSITE : Factorization Step failed with fd(x) != 1')
        return False
    if (S%2, Jacobi(Delta,N)) not in [(0,1),(1,-1)]:
        if pr: print('COMPOSITE : Jacobi step failed.')
        return False
    return True

@ trial_div([2])
def Frobenius_prob_prime(N:int, f:Polynomial, pr:bool=True):
    if not frobenius_prob_prime_core(N, f, False, pr): return False
    if pr: print('PRIME (probable) : Frobenius probable prime test passed.')
    return True

@ trial_div([2])
def Strong_Frobenius_prob_prime(N:int, f:Polynomial, pr:bool=True):
    if not frobenius_prob_prime_core(N, f, True, pr): return False
    if pr: print('PRIME (probable) : Strong Frobenius probable prime test passed.')
    return True

def quadratic_frobenius_core(N:int, b:int, c:int, pr:bool):         # core steps of quadratic frobenius test
    r1, r2 = gcd(b*b+4*c, N), gcd(c,N)
    factor = r1 if 1 < r1 < N else 0
    if 1 < r2 < N: factor = [r1,r2] if factor and r1 != r2 else r2
    if factor:
        if pr: print(HINT_FACTOR, factor)
        return False
    assert Jacobi(b*b+4*c, N) == -Jacobi(-c,N) == -1, 'Condition for Jacobi symbol not satisfied.'
    ## Step 3: compute x^{(N+1)/2} mod f(x)
    x = Polynomial(modulo=N)
    f = x*x - b*x - c
    try: y = poly_pow_mod(x, N//2+1, f)
    except ValueError as x:                               # information received which tells N is composite
        if pr: print(HINT_FACTOR, x)
        return False
    if y.degree > 0:                   # Step 3 failed
        if pr: print(f'COMPOSITE : Quadratic Frobenius test failed with (b,c) = {(b,c)} and x^((N+1)/2) not in Z/NZ.')
        return False
    ## Step 4: compute x^(N+1)
    y = y.coef[0] if y.coef else 0
    if (y*y+c) % N > 0:            # Step 4 failed
        if pr: print(f'COMPOSITE : Quadratic Frobenius test failed with (b,c) = {(b,c)} and x^(N+1) != -c (mod N,f(x)).')
        return False
    ## Step 5,6: 
    r,s = 1, N*N-1
    while s % 2 == 0:
        r += 1
        s //= 2
    try:
        y = poly_pow_mod(x,s,f)
        if y.coef == {0:1}: return True       # x^s == 1 (mod N,f(x))
        while True:
            if y.coef == {0:N-1}: return True     # x^(s*2^i) == -1 (mod N,f(x))
            r -= 1
            if r == 1: break
            y = y*y % f
    except ValueError as x:            # information received which tells N is composite
        if pr: print(HINT_FACTOR, x)
        return False
    if pr: print(f'COMPOSITE : Quadratic Frobenius test failed with (b,c) = {(b,c)} and x^s != 1 , x^(s*2^i) != -1 (mod N,f(x)).')    # Step 5 failed
    return False

@ trial_div([2,3,5,7,11,13,17,19])
@ check_sqrt
def Quadratic_Frobenius(N:int, b:int, c:int, pr:bool=True):
    if not quadratic_frobenius_core(N,b,c,pr): return False
    if pr: print('PRIME (probable) : Quadratic Frobenius test passed with (b,c) =', (b,c))
    return True

@ trial_div([2,3,5,7,11,13,17,19])
@ check_sqrt
def Quadratic_Frobenius_random(N:int, test_time:int=100, pr:bool=True):
    test_time_left = test_time
    while test_time_left > 0:                                   # randomly choose b,c
        b,c = random.randint(1,N-1), random.randint(1,N-1)
        try: ans = quadratic_frobenius_core(N,b,c,pr)    # core steps
        except AssertionError: continue
        if not ans: return False
        test_time_left -= 1
    if pr: print(f'PRIME (probable) : Quadratic Frobenius test passed with {test_time} repeated tests.')
    return True


##### 17. AKS #####

def AKS(N:int, pr:bool=True):
    assert N > 1, HINT_LARGER_THAN_1
    L = math.log2(N)
    ## Step 1
    b = 2
    while b <= L+1:
        a = powrt_int(N,b)
        if a > 0:
            if pr: print(HINT_FACTOR, a)
            return False
        b += 1
    ## Step 2,3
    r = 2
    while r < N:
        if gcd(r,N) > 1:
            if pr: print(HINT_FACTOR, gcd(r,N))
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
        if pr: print(HINT_NO_FACTOR)
        return True
    ## Step 5
    x = Polynomial(modulo=N)
    p0 = Polynomial({0:-1,r:1},N)
    xN = poly_pow_mod(x,N,p0)
    phi = Euler_phi(r)
    max_a = math.floor(math.sqrt(phi)*L)
    if pr: print(f'r = {r},  phi(r) = {phi}')
    if pr: print('max(a) =', max_a)
    if pr: print('checking if for all a, (x+a)^N == x^N + a (mod N,x^r-1) ...')
    for a in range(1,max_a+1):
        p1 = poly_pow_mod(x+a,N,p0)
        if p1.coef != (xN+a).coef:
            if pr: print(f'COMPOSITE : AKS test failed with (x+a)^N != x^N + a (mod N,x^r-1), where a,r = {a},{r}')
            return False
        #print(f'a = {a}, polynomial test passed')
    ## Step 6
    if pr: print('PRIME : AKS test passed.')
    return True
