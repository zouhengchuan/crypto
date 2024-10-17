from sage.all import *
from sympy import primerange

def is_irreducible_over_Zq(n, q):
    x = polygen(GF(q), 'x')
    polynomial = x**n - x - 1
    factors = polynomial.factor()
    return len(factors) == 1 and factors[0][1] == 1

# 定义 n 和 q 的范围
n_range = primerange(640, 800)
q_primes = list(primerange(2000, 10000))

# 存储满足条件的 (n, q) 组合
irreducible_pairs = []

# 遍历所有可能的 (n, q) 组合
for n in n_range:
    for q in q_primes:
        if is_irreducible_over_Zq(n, q):
            irreducible_pairs.append((n, q))
            print(f"Found irreducible pair: n={n}, q={q}")

# 输出所有满足条件的 (n, q) 组合
print("Irreducible pairs:")
for pair in irreducible_pairs:
    print(pair)