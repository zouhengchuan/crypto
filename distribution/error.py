from math import factorial as fac
from math import floor
import matplotlib.pyplot as plt


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
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C



p = 2048; q = 512; k = 768

fig, axs = plt.subplots(2, 2, figsize=(12, 7))

A = build_rq_law(p)
axs[0, 0].plot(A.keys(), A.values(), marker='o', markersize=1, linestyle='', color = "red")
axs[0, 0].set_title("random in R_q")

X = build_centered_binomial_law(2)
AX = law_product(A, X)

AX = iter_law_convolution(AX, k, p)
axs[0, 1].plot(AX.keys(), AX.values(), marker='o', markersize=1, linestyle='', color = "blue")
axs[0, 1].set_title("a times x")

D = {}
for i in AX:
    err = i - p/q * floor(q/p * i + 1./2)
    D[err] = D.get(err,0) + AX[i]
axs[1, 0].plot(D.keys(), D.values(), marker='o', markersize=1, linestyle='', color = "blue")
axs[1, 0].set_title("ax - p/q round(q/p ax)")

X2 = build_centered_binomial_law(2)
DX2 = law_product(X2, D)
DX2 = iter_law_convolution(DX2, k, p)
axs[1, 1].plot(DX2.keys(), DX2.values(), marker='o', markersize=1, linestyle='', color = "blue")
axs[1, 1].set_title("err * X2")

plt.suptitle('(p, q) = (%d, %d)'%(p, q), fontsize=16)
plt.tight_layout()
plt.show()