from math import factorial as fac
from math import floor

def var_of_law(D):
    div_t = 0
    avg_t = 0.
    for d in D:
        avg_t += d
    avg_t /= len(D)
    for d in D:
        div_t += D[d]*(d - avg_t)**2
    return div_t

def binomial(x, y):
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom

def centered_binomial_pdf(k, x):
    return binomial(2*k, x+k) / 2.**(2*k)

def build_centered_binomial_law(k):
    D = {}
    for i in range(-k, k+1):
        D[i] = centered_binomial_pdf(k, i)
    return D

def build_rq_law(p):
    D = {}
    for i in range(p):
        D[i] = 1./p
    return D

def law_convolution(A, B, p):
    C = {}
    for a in A:
        for b in B:
            c = (a + b) % p
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def clean_dist(A):
    B = {}
    for (x, y) in A.items():
        if y>2**(-300):
            B[x] = y
    return B

def iter_law_convolution(A, i, p):
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution(D, D, p)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A, p)
            D = clean_dist(D)
    return D

def law_product(A, B):
    C = {}
    for a in A:
        for b in B:
            c = a * b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def law_product_3n(A, B):
    C = {}
    for a in A:
        for a_prime in A:
            for b in B:
                for b_prime in B:
                    c = a * b + b_prime * (a + a_prime)
                    C[c] = C.get(c, 0) + A[a] * B[b] * A[a_prime] * B[b_prime]
    return C

def law_product_3n_2(A, B):
    C = {}
    for a in A:
        for a_prime in A:
            for b in B:
                for b_prime in B:
                    c = a * b + b_prime * a_prime
                    C[c] = C.get(c, 0) + A[a] * B[b] * A[a_prime] * B[b_prime]
    return C