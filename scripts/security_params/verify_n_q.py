import os
import re
from math import sqrt
from typing import Tuple
from MLWE_security import MLWE_summarize_attacks_test2

def process_irreducible_pairs(input_file_path: str, output_file_path: str, s):
    with open(input_file_path, 'r') as file:
        lines = file.readlines()[1:]  # 忽略第一行
    
    results = []
    for line in lines:
        match = re.match(r'\((\d+), (\d+), (\d+)\)', line.strip())
        if match:
            n, q, pq = map(int, match.groups())
            l, result = MLWE_summarize_attacks_test2(q, n, s)
            if l:
                results.append((n, q, result))
    
    with open(output_file_path, 'w') as file:
        for res in results:
            file.write(f'n={res[0]}, q={res[1]}, pq={res[2]}\n')

# 指定文件路径
s_list = []
for eta1 in range(2, 5):
    for eta2 in range(eta1, 5):
        s = sqrt(sqrt(eta1 * eta2) / 2)

        input_file_path = os.path.join('..', 'output', 'irreducible_pairs.txt')
        output_file_path = os.path.join('..', 'output', 'results_of_%d_%d.txt'% (eta1, eta2))

        # 处理不可约对并生成结果文件
        process_irreducible_pairs(input_file_path, output_file_path, s)