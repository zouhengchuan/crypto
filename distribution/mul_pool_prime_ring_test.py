from sage.all import *
from sympy import primerange
from concurrent.futures import ProcessPoolExecutor
import itertools
from MLWE_security import MLWE_summarize_attacks_test

def is_irreducible_over_Zq(n, q):
    x = polygen(GF(q), 'x')
    polynomial = x**n - x - 1
    factors = polynomial.factor()
    return len(factors) == 1 and factors[0][1] == 1

# 定义 n 和 q 的范围
# n_ranges = [primerange(640, 800), primerange(1000, 1200)]
# q_primes = list(primerange(2000, 12000))
n_ranges = [[751, 757], [761]]
q_primes = list(primerange(3000, 5000))

# 生成所有可能的 (n, q) 组合
pairs = [(n, q) for n_range in n_ranges for n in n_range for q in q_primes]

# 并行处理函数
def process_pair(pair):
    n, q = pair
    (l1, l2) = MLWE_summarize_attacks_test(q, n)
    if l1 and is_irreducible_over_Zq(n, q):
        return (n, q, l2)
    return None

# 主函数
def main():
    irreducible_pairs = []
    
    # 使用多进程池
    with ProcessPoolExecutor() as executor:
        # 使用 map 方法并行处理所有组合
        results = executor.map(process_pair, pairs)
        
        # 过滤掉 None 值
        irreducible_pairs = [result for result in results if result is not None]
        
        # 打印找到的不可约对（可选）
        for pair in irreducible_pairs:
            print(f"Found irreducible pair: n={pair[0]}, q={pair[1]}, pq={pair[2]}")
    
    # 将所有满足条件的 (n, q) 组合输出到文件
    with open('irreducible_pairs.txt', 'w') as file:
        file.write("Irreducible pairs:\n")
        for pair in irreducible_pairs:
            file.write(f"{pair}\n")
    
    # 打印完成信息
    print("All results have been written to 'irreducible_pairs.txt'.")

if __name__ == "__main__":
    main()