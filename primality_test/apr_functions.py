from .math_solve import primitive_root

def trial_division_30(N:int):
    for m in [2,3,5]:
        if m * m > N: return True
        if N % m == 0: return False
    A,k = 30,0
    I = [1,7,11,13,17,19,23,29]
    while True:
        for i in I:
            m = A*k + i
            if m == 1: continue
            if m * m > N: return True
            if N % m == 0: return False
        k += 1

# 函数v(p,q)表示q的素因数分解中p的次数
def v(p,q):
    m = 0
    while q % p == 0:
        q //= p
        m += 1
    return m

# 计算函数e(t)
def e(t):
    e = 2
    for q in range(2,t+2):
        if t % (q-1) == 0 and trial_division_30(q): e *= q**(v(q,t)+1)
    return e

# 计算q的原根并生成映射表f[x]
def generate_f(q:int):
    gx = [primitive_root(q)]
    for _ in range(q-3): gx.append(gx[0]*gx[-1] % q)
    return [0] + [gx.index(q+1-g)+1 for g in gx]