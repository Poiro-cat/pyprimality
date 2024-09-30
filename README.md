# pyprime

This package is based on Python 3, and contains prototypes of nearly all known primality test algorithms.

## List of primality test algorithms (18 in all)

### 0. Trial division

The simplest algorithm which tries to divide $N$ with all possible factors.

How to use in pyprime: `Trial_division()`, `Trial_division_30()`, `Trial_division_lim()`

### 1. Fermat's test

If there exists an integer $a$ which is coprime to $N$ and $a^{N-1}\equiv 1\pmod{N}$, then $N$ is probable prime. This is based on Fermat's little theorem.

Reference: [Fermat primality test](https://en.wikipedia.org/wiki/Fermat_primality_test)

How to use in pyprime: `Fermat()`

### 2. Miller-Rabin's test

An strong version of Fermat's test. Reference: [Miller–Rabin primality test](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test)

How to use in pyprime: `Miller_Rabin()`

### 3. Solovay–Strassen primality test

Another strong version of Fermat's test, which is related to Jacobi symbol.

Reference: [Solovay–Strassen primality test](https://en.wikipedia.org/wiki/Solovay%E2%80%93Strassen_primality_test), [Jacobi symbol](https://en.wikipedia.org/wiki/Jacobi_symbol)

How to use in pyprime: `Solovay_Strassen()`

### 4. Lucas probable prime test

A primality test based on Lucas sequence, and its strong / extra strong versions. Reference: [Lucas sequence](https://en.wikipedia.org/wiki/Lucas_sequence), [Lucas pseudoprime](https://en.wikipedia.org/wiki/Lucas_pseudoprime), [<Jon Grantham, 2000>](https://www.ams.org/journals/mcom/2001-70-234/S0025-5718-00-01197-2/S0025-5718-00-01197-2.pdf)

How to use in pyprime: `Lucas_prob_prime()`, `Strong_Lucas_prob_prime()`, `Extra_strong_Lucas_prob_prime()`

### 5. Baillie–PSW primality test

One of the most powerful probable-prime tests, which is a combination of Miller-Rabin's test based 2 and strong Lucas probable prime test.

Reference: [Baillie–PSW primality test](https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test)

How to use in pyprime: `Baillie_PSW()`

### 6. Lucas primality test

A strong version of Fermat's test using complete factorization of $N-1$. Reference: [Lucas primality test](https://en.wikipedia.org/wiki/Lucas_primality_test)

How to use in pyprime: `Lucas()`

### 7. Pocklington primality test

A weaker version of Lucas primality test, which only requires partial factorization of $N-1$.

Reference: [Pocklington primality test](https://en.wikipedia.org/wiki/Pocklington_primality_test)

How to use in pyprime: `Pocklington()`

### 8. Proth's theorem and Pépin's test (for numbers with special forms)

Proth's theorem: A primality test for Proth's numbers with forms $P=A\cdot2^n+1$, where $A$ is odd and $A< 2^n$.

Pépin's test: A primality test for Fermat numbers with forms $F_n=2^{2^n}+1$.

Reference: [Proth's theorem](https://en.wikipedia.org/wiki/Proth%27s_theorem), [Pépin's test](https://en.wikipedia.org/wiki/P%C3%A9pin%27s_test)

How to use in pyprime: `Proth()`, `Pepin()`

### 9. Factorizing Lucas probable prime test

Strong versions of Lucas probable prime test, which utilize complete or partial factorization of $N+1$.

Reference: [n+1 tests](https://t5k.org/prove/prove3_2.html)

How to use in pyprime: `Factorizing_Lucas_prob_prime_v1()`, `Factorizing_Lucas_prob_prime_v2()`

### 10. Lucas–Lehmer and Lucas–Lehmer–Riesel test (for numbers with special forms)

Lucas–Lehmer test: A primality test for Mersenne numbers with forms $M_n=2^n-1$.

Lucas–Lehmer–Riesel test: A primality test for numbers with forms $N=A\cdot2^n-1$, where $N$ is odd and $N< 2^n$.

Reference: [Lucas–Lehmer test](https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test), [Lucas–Lehmer–Riesel test](https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer%E2%80%93Riesel_test), [<Rödseth, Öystein J. 1994>](https://web.archive.org/web/20160306082833/http://folk.uib.no/nmaoy/papers/luc.pdf)

How to use in pyprime: `Lucas_Lehmer()`, `Lucas_Lehmer_Riesel()`

### 11. Goldwasser-Kilian's algorithm

The first elliptic curve prime proving (ECPP) algorithm. Reference: [Elliptic curve primality](https://en.wikipedia.org/wiki/Elliptic_curve_primality)

How to use in pyprime: `Goldwasser_Kilian()`

### 12. Atkin-Morain's algorithm

ECPP with complex multiplication method. Reference: [<A. O. L. ATKIN AND F. MORAIN, 1993>](https://www.ams.org/journals/mcom/1993-61-203/S0025-5718-1993-1199989-X/S0025-5718-1993-1199989-X.pdf)

How to use in pyprime: `Atkin_Morain()`

### 13. Tsumura's algorithm (for numbers with special forms)

ECPP for numbers with forms $N=2^k\cdot n-1$.

Reference: [<Yu Tsumura, 2009>](https://arxiv.org/pdf/0912.5279v1)

How to use in pyprime: `Tsumura()`

### 14. Methods with factors of cyclotomic polynomials

These methods are proposed by H. C. Williams et al and are related to cyclotomic polynomials of degree 3,4,6.

Reference: [<Williams, 1976, "N^2+1">](https://sci-hub.se/10.1090/s0025-5718-1976-0396390-3), [<Williams, 1976, "N^2±N+1">](https://sci-hub.se/10.1090/s0025-5718-1976-0414473-6)

How to use in pyprime: `Quartic_field()`, `Cubic_Sextic_field()`

### 15. APR-CL algorithm

The generalization of H. C. Williams's methods. This is the first **deterministic** primality test within **nearly polynomial time** which is **applicable for any integers**.

Reference: [<APR, 1983>](https://www.jstor.org/stable/2006975), [<APR-CL, 1984>](https://www.ams.org/journals/mcom/1984-42-165/S0025-5718-1984-0726006-X/S0025-5718-1984-0726006-X.pdf)

How to use in pyprime: `APR_CL()`

### 16. Frobenius probable prime test

A series of primality tests using polynomials over finite fields. Reference: [<Jon Grantham, 2000>](https://www.ams.org/journals/mcom/2001-70-234/S0025-5718-00-01197-2/S0025-5718-00-01197-2.pdf), [Quadratic Frobenius test](https://en.wikipedia.org/wiki/Quadratic_Frobenius_test)

How to use in pyprime: `Frobenius_prob_prime()`, `Strong_Frobenius_prob_prime()`, `Quadratic_Frobenius()`, `Random_Quadratic_Frobenius()`

### 17. AKS test

Based on the Fermat's little theorem over polynomial rings, AKS test is the first **deterministic** primality test which is able to determine **in polynomial time**, whether an **arbitrary given number** is prime or composite and this **without relying on mathematical conjectures** (such as generalized Riemann hypothesis that the deterministic version of Miller's test relies on).

It's very efficient according to its complexity, theoretically, but is extremely slow in practice.

Reference: [<AKS, 2004>](https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf), [AKS primality test](https://en.wikipedia.org/wiki/AKS_primality_test)

How to use in pyprime: `AKS()`

## Other algorithms

`pyprime` also supports other algorithms, such as the following:

### Operation between polynomials

One could define a series operations between polynomials (over finite fields), including addition, subtraction, multiplication, modular division, greatest common divisor, et al.

How to use in pyprime: `Polynomial()`

### Group operation in elliptic curve group over finite fields

An elliptic curve over a finite field forms an Abelian group, where addition operation is defined between points. One can also define multiplication operation between a point and a natural number.

Reference: [Elliptic curve](https://en.wikipedia.org/wiki/Elliptic_curve)

How to use in pyprime: `Elliptic_curve().plus()`, `Elliptic_curve().times()`

### Schoof's algorithm

An algorithm which computes the order of an elliptic curve group over finite fields. Reference: [Schoof's algorithm](https://en.wikipedia.org/wiki/Schoof%27s_algorithm)

How to use in pyprime: `Elliptic_curve().schoof()`
