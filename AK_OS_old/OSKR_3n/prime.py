def is_prime(n):
    """Check if the given number is prime."""
    if n <= 1:
        return False
    if n <= 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

divisor = 3**1 * 2**6  # 3^2 * 2^7 = 9 * 128 = 1152

primes = []
for q in range(1000, 20001):  # Check numbers from 1000 to 4000 inclusive
    if is_prime(q) and (q - 1) % divisor == 0:
        primes.append(q)

print("Primes q in the range [1000, 4000] where q-1 is divisible by 3^2 * 2^7:")
print(primes)

# 3^2 * 2^6 [1153, 3457, 6337, 7489, 8641, 10369, 12097, 13249, 14401, 18433, 19009]

# 3^1 * 2^7 [1153, 2689, 3457, 4993, 6529, 7297, 7681, 9601, 10369, 10753, 12289, 13441, 14593, 15361, 18049, 18433]

# 3^1 * 2^8 [7681, 10753, 12289, 14593, 15361, 18433]