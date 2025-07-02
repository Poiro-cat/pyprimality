# pyprimality

This is a Python 3 package implementing prototypes of nearly all known primality test algorithms.

## List of primality test algorithms (18 types in all)

### 0. Trial division

The simplest algorithm which checks divisibility by all possible factors.

- `Trial_division(N)` : Trial divisors range from 2 to $\lfloor\sqrt N\rfloor$
- `Trial_division_30(N)` : After testing divisibility by $2,3,5$, trial devisors are chosen  by $m=30k+i$, where $k\geq0$ and $i$ is coprime with $30$.
- `Trial_division_lim(N)` : Trial devisors are chosen to from a dynamically generated prime number list.

### 1. Fermat's test

If there exists an integer $a$ which is coprime to $N$ and $a^{N-1}\equiv 1\pmod{N}$, then $N$ is a probable prime. This is based on Fermat's little theorem.

Reference: [\<Fermat primality test\>](https://en.wikipedia.org/wiki/Fermat_primality_test)

- `Fermat(N)` : Fermat's test for $N$ with a set of random chosen base numbers.

### 2. Miller-Rabin's test

A stronger version of Fermat's test with additional verification conditions.

Reference: [\<Miller-Rabin primality test\>](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test)

- `Miller_Rabin(N)` : Miller-Rabin's test for $N$ with a set of random chosen base numbers.

### 3. Solovay-Strassen primality test

A stronger version of Fermat's test, which is related to Jacobi symbol.

Reference: [\<Solovay-Strassen primality test\>](https://en.wikipedia.org/wiki/Solovay%E2%80%93Strassen_primality_test), [\<Jacobi symbol\>](https://en.wikipedia.org/wiki/Jacobi_symbol)

- `Solovay_Strassen(N)` : Solovay-Strassen's test for $N$ with a set of random chosen base numbers.

### 4. Lucas probable prime test

A type of primality test based on Lucas sequences. If certain identities about Lucas sequences hold for the tested integer $N$, then $N$ is a probable prime. Several different versions of Lucas probable prime test are presented in `pyprimality`.

Reference: [\<Lucas sequence\>](https://en.wikipedia.org/wiki/Lucas_sequence), [\<Lucas pseudoprime\>](https://en.wikipedia.org/wiki/Lucas_pseudoprime), [\<Jon Grantham, 2000\>](https://www.ams.org/journals/mcom/2001-70-234/S0025-5718-00-01197-2/S0025-5718-00-01197-2.pdf)

- `Lucas_prob_prime(N, P, Q)` : Lucas probable prime test for $N$ using $P,Q$ as parameters.
- `Strong_Lucas_prob_prime(N, P, Q)` : Strong version of the test for $N$ using $P,Q$ as parameters.
- `Extra_strong_Lucas_prob_prime(N, b)` : Extra-strong version of the test for $N$ using $b$ as the parameter.
- `Lucas_prob_prime_random(N)` : Lucas probable prime test for $N$ using a set of randomly chosen parameters.
- `Strong_Lucas_prob_prime_random(N)` : Strong version of the test for $N$ using a set of randomly chosen parameters.
- `Extra_strong_Lucas_prob_prime_random(N)` : Extra-strong version of the test for $N$ using a set of randomly chosen parameters.

### 5. Baillie-PSW primality test

One of the most powerful probable-prime tests, a combination of Miller-Rabin's test (base 2) and strong Lucas probable prime test (with specially chosen parameters).

Reference: [\<Baillie-PSW primality test\>](https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test)

- `Baillie_PSW(N)` : Baillie-PSW test for $N$.

### 6. Lucas primality test

A stronger version of Fermat's test using complete factorization of $N-1$. Reference: [\<Lucas primality test\>](https://en.wikipedia.org/wiki/Lucas_primality_test)

- `Lucas(N)` : Lucas primality test for $N$, which automatically factorizes $N-1$ by trial division.

### 7. Pocklington primality test

A variant of Lucas primality test, which only requires partial factorization of $N-1$.

Reference: [\<Pocklington primality test\>](https://en.wikipedia.org/wiki/Pocklington_primality_test)

- `Pocklington(N, known_factors)` : Pocklington test for $N$. You may optionally provide known prime factors of $N-1$ as `known_factors` to assist the factorization.

### 8. Proth's theorem and Pépin's test (for numbers with special forms)

Proth's theorem: A primality test for Proth's numbers with forms $P=A\cdot2^n+1$ where $A$ is odd and $A< 2^n$.

Pépin's test: A primality test for Fermat numbers with forms $F_n=2^{2^n}+1$.

Reference: [\<Proth's theorem\>](https://en.wikipedia.org/wiki/Proth%27s_theorem), [\<Pépin's test\>](https://en.wikipedia.org/wiki/P%C3%A9pin%27s_test)

- `Proth(A, n)` : Proth's test for $A\cdot2^n+1$.
- `Pepin(n)` : Pépin's test for $2^{2^n}+1$.

### 9. Factorizing Lucas probable prime test

Strong versions of Lucas probable prime test, which utilize complete or partial factorization of $N+1$.

Reference: [\<n+1 tests\>](https://t5k.org/prove/prove3_2.html)

- `Factorizing_Lucas_prob_prime_v0(N)` : Factorizing Lucas probable prime test for $N$, using complete factorization of $N+1$.
- `Factorizing_Lucas_prob_prime_v1(N, known_factors)` : Factorizing Lucas probable prime test for $N$, using partial factorization of $N+1$. You may optionally provide known prime factors of $N+1$ as `known_factors` to assist the factorization.

### 10. Lucas-Lehmer and Lucas-Lehmer-Riesel test (for numbers with special forms)

Lucas-Lehmer test: A primality test for Mersenne numbers with forms $M_n=2^n-1$.

Lucas-Lehmer-Riesel test: A primality test for numbers with forms $N=A\cdot2^n-1$ where $A$ is odd and $A< 2^n$.

Reference: [\<Lucas-Lehmer test\>](https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer_primality_test), [\<Lucas-Lehmer-Riesel test\>](https://en.wikipedia.org/wiki/Lucas%E2%80%93Lehmer%E2%80%93Riesel_test), [\<Rödseth, Öystein J. 1994\>](https://web.archive.org/web/20160306082833/http://folk.uib.no/nmaoy/papers/luc.pdf)

- `Lucas_Lehmer(p)` : Lucas-Lehmer test for $2^n-1$.
- `Lucas_Lehmer_Riesel(A, n)` : Lucas-Lehmer-Riesel test for $A\cdot2^n-1$.

### 11. Goldwasser-Kilian's algorithm

The first elliptic curve primality proving (ECPP) algorithm. ECPP is a type of primality test which utilizes the group order of elliptic curves over a finite field.

Reference: [\<Elliptic curve primality\>](https://en.wikipedia.org/wiki/Elliptic_curve_primality)

- `Goldwasser_Kilian(N)` : The algorithm generates elliptic curves by randomly choosing parameters.

### 12. Atkin-Morain's algorithm

ECPP with complex multiplication method. Reference: [\<A. O. L. ATKIN AND F. MORAIN, 1993\>](https://www.ams.org/journals/mcom/1993-61-203/S0025-5718-1993-1199989-X/S0025-5718-1993-1199989-X.pdf)

- `Atkin_Morain(N)` : The algorithm generates elliptic curves with certain group orders by complex multiplication method.

### 13. Tsumura's algorithm (for numbers with special forms)

ECPP for numbers with forms $N=2^k\cdot n-1$.

Reference: [\<Yu Tsumura, 2009\>](https://arxiv.org/pdf/0912.5279v1)

- `Tsumura(k, n)` : Tsumura's ECPP for $2^k\cdot n-1$.

### 14. Methods with factors of cyclotomic polynomials

These methods are proposed by H. C. Williams et al and utilize factorization of cyclotomic polynomials of degree 3,4,6.

Reference: [\<Williams, 1976, "N^2+1"\>](https://sci-hub.se/10.1090/s0025-5718-1976-0396390-3), [\<Williams, 1976, "N^2±N+1"\>](https://sci-hub.se/10.1090/s0025-5718-1976-0414473-6)

- `Cubic_field(N)` : Primality test for $N$ using partial factorization of $N^2+N+1$.
- `Quartic_field(N)` : Primality test for $N$ using partial factorization of $N^2+1$.
- `Sextic_field(N)` : Primality test for $N$ using partial factorization of $N^2-N+1$.

### 15. APR-CL algorithm

The generalization of H. C. Williams's methods. This is the first **deterministic** primality test within **nearly polynomial time** which is **applicable for any integers**.

Reference: [\<APR, 1983\>](https://www.jstor.org/stable/2006975), [\<APR-CL, 1984\>](https://www.ams.org/journals/mcom/1984-42-165/S0025-5718-1984-0726006-X/S0025-5718-1984-0726006-X.pdf)

- `APR_CL(N)` : APR-CL test for $N$.

### 16. Frobenius probable prime test

A type of primality test by checking several properties of a polynomial over a finite field.

Reference: [\<Jon Grantham, 2000\>](https://www.ams.org/journals/mcom/2001-70-234/S0025-5718-00-01197-2/S0025-5718-00-01197-2.pdf), [\<Quadratic Frobenius test\>](https://en.wikipedia.org/wiki/Quadratic_Frobenius_test)

- `Frobenius_prob_prime(N, f)` : Frobenius probable prime test for $N$ with polynomial $f(x)$. The type of input `f` is a class defined in `pyprime` (see "Polynomials" in this file).
- `Strong_Frobenius_prob_prime(N, f)` : Stronger version of Frobenius probable prime test for $N$ with polynomial $f(x)$.
- `Quadratic_Frobenius(N, b, c)` : Special version of Frobenius probable prime test for $N$ with quadratic polynomial $f(x)=x^2-bx-c$.
- `Quadratic_Frobenius_random(N)` : Quadratic Frobenius test for $N$ with randomly chosen $b,c$.

### 17. AKS test

Based on the Fermat's little theorem over polynomial rings, AKS test is the first **deterministic** primality test which is able to determine **in polynomial time**, whether an **arbitrary given number** is prime or composite and this **without relying on mathematical conjectures** (such as generalized Riemann hypothesis that the deterministic version of Miller's test relies on).

It's very efficient according to its complexity, theoretically, but is extremely slow in practice.

Reference: [\<AKS, 2004\>](https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf), [\<AKS primality test\>](https://en.wikipedia.org/wiki/AKS_primality_test)

- `AKS(N)` : AKS test for $N$.

## Other

`pyprimality` also provides other algorithms, such as the following:

### Polynomials

One could define a series operations between polynomials (over finite fields), including addition, subtraction, multiplication, modular division, greatest common divisor, et al.

- `Polynomial(coef, modulo)` : It defines a polynomial with given coefficients. The type of input `coef` should be `list` or `dict`. You may optionally provide a positive number as `modulo`, which means the polynomial is defined over a finite field whose characteristic is the given number.

### Elliptic curve group over finite fields

An elliptic curve over a finite field forms an Abelian group, with addition defined for points and scalar multiplication defined for points and natural numbers.

Reference: [\<Elliptic curve\>](https://en.wikipedia.org/wiki/Elliptic_curve)

- `Point(x, y)` : It defines an integer point $(x, y)$ in two-dimensional Cartesian coordinate. Notice: If $x=y=0$, then it is the point at infinity instead of the origin.
- `Elliptic_curve(N, a, b)` : It defines an elliptic curve over $\mathbb Z/N\mathbb Z$ with equation $y^2=x^3+ax+b\pmod N$. 
- `Elliptic_curve(N, a, b).plus(P, Q)` : It returns the result of group operation $P+Q$, where $P,Q$ are both points on the elliptic curve.
- `Elliptic_curve(N, a, b).times(m, P)` : It returns the result of scalar multiplication $mP$, where $m\in\mathbb N$ and $P$ is a point on the elliptic curve.

### Schoof's algorithm

An algorithm which computes the order of an elliptic curve group over finite fields. Reference: [\<Schoof's algorithm\>](https://en.wikipedia.org/wiki/Schoof%27s_algorithm)

- `Elliptic_curve(N, a, b).schoof()` : It returns the group order of the elliptic curve.
