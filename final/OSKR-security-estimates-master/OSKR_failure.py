import operator as op
from math import factorial as fac
from math import sqrt, log
import sys
from proba_util import *
from scipy.stats import norm

# 3n 上的乘法
def error_probability(ps):
    F = error_distribution(ps)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    f = ps.n * proba
    # f = 0
    # for k in range(ps.n):
    #     F = error_distribution(ps, k)
    #     proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    #     f += proba

    print("mul in 3n")
    print("failure origin: = 2^%.5f"% (log(f + 2.**(-300))/log(2)))

def error_distribution(ps, k = None):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    if k is None:
        k = ps.n // 2
    chisr = build_centered_binomial_law(ps.ks)                # LWE error law for the key (s, r)
    chie = build_centered_binomial_law(ps.ke)                 # LWE error law for the error (e)
    chie_ct = build_centered_binomial_law(ps.ke_ct)           # LWE error law for the ciphertext (e1,e2)

    Ru = build_mod_switching_error_law(ps.q, ps.rqc)          # rounding error first ciphertext (eu)
    chiRe = law_convolution_q(chie_ct, Ru, ps.q)                      # LWE + rounding error ciphertext  (e1 + eu)

    if (k < ps.n // 2):
        B1 = iter_law_product_3n(chie, chisr, ps.n, ps.q, k+1)  # (LWE+Rounding error) * LWE (as in a E*S product) (e*r)
        B2 = iter_law_product_3n(chisr, chiRe, ps.n, ps.q, k+1) # (s*e1 + s*eu)
    else:
        B1 = iter_law_product_3n(chie, chisr, ps.n, ps.q, 0)
        B2 = iter_law_product_3n(chisr, chiRe, ps.n, ps.q, 0)

    C1 = iter_law_convolution_q(B1, ps.m, ps.q)                         # m*(e*r)
    C2 = iter_law_convolution_q(B2, ps.m, ps.q)                         # m*(s*e1 + s*eu) 

    C = law_convolution_q(C1, C2, ps.q)                               # m*n*(e*r) + m*n*(s*e1+s*eu)

    Rv = build_mod_switching_error_law(ps.q, ps.rq2)          # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution_q(Rv, chie, ps.q)                             # LWE+Rounding2 error (ev + e2)

    D = law_convolution_q(C, F, ps.q)                                 # Final error (m*n*(e*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D


# 直接按项数估计
def error_probability_1(ps):
    F = error_distribution_1(ps)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    f = ps.n * proba
    # f = 0
    # for k in range(ps.n):
    #     F = error_distribution_1(ps, k)
    #     proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    #     f += proba

    print("3n/2 times convolution")
    print("failure origin: = 2^%.5f"% (log(f + 2.**(-300))/log(2)))

def error_distribution_1(ps, k = None):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    if k is None:
        k = ps.n // 2

    if k < ps.n // 2:
        n = ps.n * 3 // 2 - k - 1
    else:
        n = ps.n * 3 // 2

    chisr = build_centered_binomial_law(ps.ks)          # LWE error law for the key (s, r)
    chie = build_centered_binomial_law(ps.ke)           # LWE error law for the error (e)
    chie_ct = build_centered_binomial_law(ps.ke_ct)     # LWE error law for the ciphertext (e1,e2)

    Ru = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRe = law_convolution(chie_ct, Ru)                # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie, chisr)                       # (LWE+Rounding error) * LWE (as in a E*S product) (e*r)
    B2 = law_product(chisr, chiRe)                      # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, round(ps.m * n))      # m*n*(e*r)
    C2 = iter_law_convolution(B2, round(ps.m * n))      # m*n*(s*e1 + s*eu) 

    C = law_convolution(C1, C2)                         # m*n*(e*r) + m*n*(s*e1+s*eu)

    Rv = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(Rv, chie)                       # LWE+Rounding2 error (ev + e2)

    D = law_convolution(C, F)                           # Final error (m*n*(e*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D


# 展开到二次(in 3n)
def error_probability_2(ps):
    proba1 = {}
    proba2 = {}
    D = error_distribution(ps, 0)
    proba1[0] = tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
    proba2[0] = 0

    for k in range(1, ps.n//2 + 1):
        D = error_distribution(ps, k)
        d = tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
        proba1[k] = proba1[k-1] + d
        proba2[k] = proba2[k-1] + proba1[k-1] * d
    for k in range(ps.n//2 + 1, ps.n):
        proba1[k] = proba1[k-1] + d
        proba2[k] = proba2[k-1] + proba1[k-1] * d

    f = proba1[ps.n-1] - proba2[ps.n-1]
    print(proba1[ps.n-1])
    print(proba2[ps.n-1])

    print("mul in 3n of 2")
    print("failure origin: = 2^%.5f"% (log(f + 2.**(-300))/log(2)))

# 展开到二次(按项数)
def error_probability_3(ps):
    E, F = error_distribution_3(ps)

    proba1 = {}
    proba2 = {}
    C = iter_law_convolution(E, ps.m * ps.n)
    D = law_convolution(C,F)
    proba1[0] = tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
    proba2[0] = 0

    for k in range(1,ps.n//2+1):
        C = iter_law_convolution(E, ps.m * (k + ps.n))
        D = law_convolution(C,F)
        d = tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
        proba1[k] = proba1[k-1] + d
        proba2[k] = proba2[k-1] + proba1[k-1] * d
    for k in range(ps.n//2+1,ps.n):
        proba1[k] = proba1[k-1] + d
        proba2[k] = proba2[k-1] + proba1[k-1] * d

    f = proba1[ps.n-1] - proba2[ps.n-1]
    print(proba1[ps.n-1])
    print(proba2[ps.n-1])

    print("3n/2 times convolution of 2")
    print("failure origin: = 2^%.5f"% (log(f + 2.**(-300))/log(2)))

def error_distribution_3(ps):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chisr = build_centered_binomial_law(ps.ks)          # LWE error law for the key (s, r)
    chie_ct = build_centered_binomial_law(ps.ke_ct)     # LWE error law for the ciphertext (e1,e2)
    chie = build_centered_binomial_law(ps.ke)           # (e)

    Ru = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRe = law_convolution(chie_ct, Ru)                # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie, chisr)                       # (LWE+Rounding error) * LWE (as in a E*S product) (e*r)
    B2 = law_product(chisr, chiRe)                      # (s*e1 + s*eu)
   
    E = law_convolution(B1, B2)

    Rv = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(Rv, chie)                       # LWE+Rounding2 error (ev + e2)

    return E, F


def bandwidth_K320(ps):
    seedLen = 40
    pk = seedLen + ceil((ps.m * ps.n * ceil(log(ps.q, 2))) / 8)
    ct = ceil(ps.m * ps.n * ceil(log(ps.rqc, 2)) / 8) + ceil(320 * ceil(log(ps.rq2, 2)) / 8)
    print('|K| = %d, |pk| = %d, |ct| = %d, bandwidth = %d\n'%(seedLen, pk, ct, (pk+ct)))


def bandwidth_K512(ps):
    seedLen = 64
    pk = seedLen + ceil((ps.m * ps.n * ceil(log(ps.q, 2))) / 8)
    ct = ceil(ps.m * ps.n * ceil(log(ps.rqc, 2)) / 8) + ceil(512 * ceil(log(ps.rq2, 2)) / 8)
    print('|K| = %d, |pk| = %d, |ct| = %d, bandwidth = %d\n'%(seedLen, pk, ct, (pk+ct)))