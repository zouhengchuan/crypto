import operator as op
from math import factorial as fac
from math import sqrt, log
import sys
from tqdm import tqdm
import time
from proba_util import *

def p2_cyclotomic_final_error_distribution(ps):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
    chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
    chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)
    Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)
    Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key  (e + et)
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r + et*r)
    B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, 2 * ps.m * ps.n)          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, 2 * ps.m * ps.n)          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D


def p2_cyclotomic_final_error_distribution_AKCN(ps):    # Decompress rounding -> floor
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
    chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
    chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)
    Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)
    Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key  (e + et)
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r + et*r)
    B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, 2 * ps.m * ps.n)          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, 2 * ps.m * ps.n)          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law_AKCN(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D

def p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps,i):    # Decompress rounding -> floor
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
    chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
    chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)
    Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)
    Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key  (e + et)
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r + et*r)
    B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, ps.m * i)          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, ps.m * i)          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law_AKCN(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D

def p2_cyclotomic_final_error_distribution_floor(ps, i):    # Decompress rounding -> floor
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
    chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
    chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)
    Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)
    Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key  (e + et)
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r + et*r)
    B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, 2 * ps.m * ps.n)          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, 2 * ps.m * ps.n)          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law_floor(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D


def p2_cyclotomic_error_probability(ps):
    F = p2_cyclotomic_final_error_distribution(ps)
    #proba = tail_probability(F, ps.q/4)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    return F, ps.n*proba

def p2_cyclotomic_error_probability_AKCN(ps):      # Decompress rounding -> AKCN
    F = p2_cyclotomic_final_error_distribution_AKCN(ps)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    return F, ps.n*proba

def p2_cyclotomic_error_probability_AKCN_onebyone_all(ps):      # Decompress rounding -> AKCN all degree
    proba = {}
    F = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps,ps.n)
    proba[0] = 1. - tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    for i in range(1, ps.n-1):
        D = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps, 2*ps.n-i)
        proba[i] = proba[i-1] * (1 - tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5))
    G = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps, ps.n+1)
    proba[ps.n-1] = proba[ps.n-2] * (1 - tail_probability_frac(G, -ps.q/4., ps.q/4. - 0.5))
    return (1 - proba[ps.n-1])

def p2_cyclotomic_error_probability_AKCN_onebyone_1(ps):      # Decompress rounding -> AKCN 1 degree
    proba = {}
    F = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps,ps.n)
    proba[0] = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    for i in range(1, ps.n-1):
        D = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps, 2*ps.n-i)
        proba[i] = proba[i-1] + tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
    G = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps, ps.n+1)
    proba[ps.n-1] = proba[ps.n-2] + tail_probability_frac(G, -ps.q/4., ps.q/4. - 0.5)
    return proba[ps.n-1]

def p2_cyclotomic_error_probability_AKCN_onebyone_2(ps):      # Decompress rounding -> AKCN 2 degree
    proba1 = {}
    proba2 = {}
    F = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps,ps.n)
    proba1[0] = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    proba2[0] = 0
    for i in tqdm(range(1, ps.n-1)):
        D = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps, 2*ps.n-i)
        proba1[i] = proba1[i-1] + tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
        proba2[i] = proba2[i-1] + proba1[i-1] * tail_probability_frac(D, -ps.q/4., ps.q/4. - 0.5)
    G = p2_cyclotomic_final_error_distribution_AKCN_onebyone(ps, ps.n+1)
    proba1[ps.n-1] = proba1[ps.n-2] + tail_probability_frac(G, -ps.q/4., ps.q/4. - 0.5)
    proba2[ps.n-1] = proba2[ps.n-2] + proba1[ps.n-2] * tail_probability_frac(G, -ps.q/4., ps.q/4. - 0.5)
    return (proba1[ps.n-1] - proba2[ps.n-1])

def p2_cyclotomic_error_probability_floor(ps):      # Decompress rounding -> floor
    F = p2_cyclotomic_final_error_distribution_floor(ps)
    #proba = tail_probability(F, ps.q/4.)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    return F, ps.n*proba