from math import factorial as fac
from math import floor
import matplotlib.pyplot as plt
from prob import *


def var_p(AX, q, p, str):
    D = {}
    for i in AX:
        err = i - q/p * floor(p/q * i + 1./2)
        D[err] = D.get(err,0) + AX[i]
    var = var_of_law(D)
    print("(prime) var of %s = %.20f"% (str, var))



q = 1024; p = 512; k = 768

D = {}
for i in range(q):
    err = i - q/p * floor(p/q * i + 1./2)
    D[err] = D.get(err,0) + 1./q
var = var_of_law(D)
print("var of uni = %.20f"% var)


A = build_rq_law(q)
X = build_centered_binomial_law(2)
AX_1 = law_product(A, X)

AX = iter_law_convolution(AX_1, k, q)
var_p(AX, q, p, "0")

for j in range(1, k):
    AX = iter_law_convolution(AX_1, (2 * k - j), q)
    var_p(AX, q, p, "%d"% j)



# fig, axs = plt.subplots(2, 2, figsize=(12, 7))
# axs[0, 0].plot(A.keys(), A.values(), marker='o', markersize=1, linestyle='', color = "red")
# axs[0, 0].set_title("random in R_q")
# axs[0, 1].plot(AX.keys(), AX.values(), marker='o', markersize=1, linestyle='', color = "blue")
# axs[0, 1].set_title("a times x")
# axs[1, 0].plot(D.keys(), D.values(), marker='o', markersize=1, linestyle='', color = "blue")
# axs[1, 0].set_title("p/q*ax - round(p/q ax)")
# axs[1, 1].plot(DX2.keys(), DX2.values(), marker='o', markersize=1, linestyle='', color = "blue")
# axs[1, 1].set_title("err * X2")

# plt.suptitle('(p, q) = (%d, %d)'%(p, q), fontsize=16)
# plt.tight_layout()
# plt.show()