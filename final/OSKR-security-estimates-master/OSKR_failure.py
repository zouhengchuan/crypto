import operator as op
from math import factorial as fac
from math import sqrt, log
import sys
from proba_util import *
from scipy.stats import norm

def p2_cyclotomic_error_probability(ps):
    F = p2_cyclotomic_final_error_distribution(ps)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    f = ps.n * proba

    print("test of 3n")
    print("failure origin: %.1f = 2^%.5f"%(f, log(f + 2.**(-300))/log(2)))

def p2_cyclotomic_final_error_distribution(ps,k = None):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    if k is None:
        k = ps.n * 3 // 2
    chisr = build_centered_binomial_law(ps.ks)          # LWE error law for the key (s, r)
    chie = build_centered_binomial_law(ps.ke)           # LWE error law for the error (e)
    chie_ct = build_centered_binomial_law(ps.ke_ct)     # LWE error law for the ciphertext (e1,e2)

    Ru = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRe = law_convolution(chie_ct, Ru)                # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie, chisr)                       # (LWE+Rounding error) * LWE (as in a E*S product) (e*r)
    B2 = law_product(chisr, chiRe)                      # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, round(ps.m * k))      # m*n*(e*r)
    C2 = iter_law_convolution(B2, round(ps.m * k))      # m*n*(s*e1 + s*eu) 

    C = law_convolution(C1, C2)                         # m*n*(e*r) + m*n*(s*e1+s*eu)

    Rv = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(Rv, chie)                       # LWE+Rounding2 error (ev + e2)

    D = law_convolution(C, F)                           # Final error (m*n*(e*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D


# 展开到二次
def p2_cyclotomic_error_probability_2(ps):
    E, F = p2_cyclotomic_final_error_distribution_2(ps)

    proba1 = {}
    proba2 = {}
    C = iter_law_convolution(E, ps.m * ps.n)
    D = law_convolution(C,F)
    proba1[0] = tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
    proba2[0] = 0

    for k in range(1,ps.n//2+1):
        C = iter_law_convolution(E, ps.m * (k + ps.n))
        D = law_convolution(C,F)
        proba1[k] = proba1[k-1] + tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
        proba2[k] = proba2[k-1] + proba1[k-1] * tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
    for k in range(ps.n//2+1,ps.n):
        proba1[k] = proba1[k-1] + tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
        proba2[k] = proba2[k-1] + proba1[k-1] * tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
    proba = proba1[ps.n-1] - proba2[ps.n-1]
    print(proba1[ps.n-1])
    print(proba2[ps.n-1])

    f = tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5) * ps.n
    print ("failure origin: %.1f = 2^%.5f"%(f, log(f + 2.**(-300))/log(2)))

    print("test2 of 3n")
    return F, proba

def p2_cyclotomic_final_error_distribution_2(ps):
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
    ct = ceil(ps.m * ps.n * ceil(log(ps.rqc, 2)) / 8) + ceil(ps.n * ceil(log(ps.rq2, 2)) / 8)
    print('|K| = %d, |pk| = %d, |ct| = %d, bandwidth = %d\n'%(seedLen, pk, ct, (pk+ct)))


def bandwidth_K512(ps):
    seedLen = 64
    pk = seedLen + ceil((ps.m * ps.n * ceil(log(ps.q, 2))) / 8)
    ct = ceil(ps.m * ps.n * ceil(log(ps.rqc, 2)) / 8) + ceil(ps.n * ceil(log(ps.rq2, 2)) / 8)
    print('|K| = %d, |pk| = %d, |ct| = %d, bandwidth = %d\n'%(seedLen, pk, ct, (pk+ct)))