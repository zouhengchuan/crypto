import operator as op
from math import factorial as fac
from math import sqrt, log
import sys
from proba_util import *
from scipy.stats import norm

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

    C1 = iter_law_convolution(B1, ps.m * ps.n)          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, ps.m * ps.n)          # m*n*(s*e1 + s*eu) 

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

    C1 = iter_law_convolution(B1, ps.m * ps.n)          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, ps.m * ps.n)          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law_AKCN(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D



def p2_cyclotomic_error_probability(ps):
    F = p2_cyclotomic_final_error_distribution_test(ps)
    #proba = tail_probability(F, ps.q/4)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    return F, ps.n*proba

def p2_cyclotomic_error_probability_AKCN(ps):      # Decompress rounding -> floor
    F = p2_cyclotomic_final_error_distribution_AKCN(ps)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    return F, ps.n*proba


def p2_cyclotomic_final_error_distribution_test(ps):
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

    C1 = iter_law_convolution(B1, round(ps.m * ps.n * 1.5))          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, round(ps.m * ps.n * 1.5))          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
    print("test of 3n")
    return D


def error_rate_test(ps):
    div_t = 0
    avg_t = 0
    for i in range(ps.q):
        temp = floor(ps.q / ps.rqc * floor(ps.rqc / ps.q * i + 1/2) + 1/2) - i
        avg_t += temp / ps.q
    for i in range(ps.q):
        temp = floor(ps.q / ps.rqc * floor(ps.rqc / ps.q * i + 1/2) + 1/2) - i
        div_t += (temp-avg_t)**2 / ps.q

    s0 = 1.5*ps.n*ps.m*ps.ks/2*(ps.ke+div_t)+ps.ke/2
    dis = ps.q/4-1-ps.q/(2*2**ps.rq2)
    s = sqrt(s0)
    pr = norm.logsf(dis/s) + 1
    print("err of 3n:")
    print("    = 2^%.2f"% pr)
