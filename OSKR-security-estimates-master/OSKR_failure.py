import operator as op
from math import factorial as fac
from math import sqrt, log
import sys
from proba_util import *
from scipy.stats import norm

# def p2_cyclotomic_final_error_distribution(ps):
#     """ construct the final error distribution in our encryption scheme
#     :param ps: parameter set (ParameterSet)
#     """
#     chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
#     chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
#     chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)
#     Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)
#     Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
#     chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key  (e + et)
#     chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

#     B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r + et*r)
#     B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

#     C1 = iter_law_convolution(B1, ps.m * ps.n)          # m*n*(e*r + et*r) prime*2
#     C2 = iter_law_convolution(B2, ps.m * ps.n)          # m*n*(s*e1 + s*eu) 

#     C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

#     R2 = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
#     F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
#     D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
#     return D

# def p2_cyclotomic_final_error_distribution_AKCN(ps):    # Decompress rounding -> floor
#     """ construct the final error distribution in our encryption scheme
#     :param ps: parameter set (ParameterSet)
#     """
#     chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
#     chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
#     chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)
#     Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)
#     Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
#     chiRs = law_convolution(chis, Rk)                   # LWE+Rounding error key  (e + et)
#     chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

#     B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r + et*r)
#     B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

#     C1 = iter_law_convolution(B1, ps.m * ps.n)          # m*n*(e*r + et*r) prime*2
#     C2 = iter_law_convolution(B2, ps.m * ps.n)          # m*n*(s*e1 + s*eu) 

#     C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

#     R2 = build_mod_switching_error_law_AKCN(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
#     F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
#     D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
#     return D


# def p2_cyclotomic_error_probability(ps):
#     F = p2_cyclotomic_final_error_distribution(ps)
#     #proba = tail_probability(F, ps.q/4)
#     proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
#     return F, ps.n*proba

# def p2_cyclotomic_error_probability_AKCN(ps):      # Decompress rounding -> floor
#     F = p2_cyclotomic_final_error_distribution_AKCN(ps)
#     proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
#     return F, ps.n*proba


def p2_cyclotomic_error_probability_test(ps):
    F = p2_cyclotomic_final_error_distribution_test(ps)
    proba = tail_probability_frac(F, -ps.q/4., ps.q/4. - 0.5)
    print("test of 3n")
    return F, ps.n*proba

def p2_cyclotomic_final_error_distribution_test(ps,k = None):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    if k is None:
        k = ps.n * 3 // 2
    chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
    chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
    chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)

    # Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)

    Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRs = chie_pk                                     # LWE+Rounding error key  (e)
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r)
    B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, round(ps.m * k))          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, round(ps.m * k))          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)
    return D


def p2_cyclotomic_error_probability_test2(ps):
    E, F = p2_cyclotomic_final_error_distribution_test2(ps)

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

def p2_cyclotomic_final_error_distribution_test2(ps):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chis = build_centered_binomial_law(ps.ks)           # LWE error law for the key (s)
    chie = build_centered_binomial_law(ps.ke_ct)        # LWE error law for the ciphertext (e1,e2)
    chie_pk = build_centered_binomial_law(ps.ke)        # (e,r)

    # Rk = build_mod_switching_error_law(ps.q, ps.rqk)    # Rounding error public key (et)

    Rc = build_mod_switching_error_law(ps.q, ps.rqc)    # rounding error first ciphertext (eu)
    chiRs = chie_pk                                     # LWE+Rounding error key  (e)
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chie_pk, chiRs)                    # (LWE+Rounding error) * LWE (as in a E*S product) (e*r)
    B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)
   
    E = law_convolution(B1, B2)

    # C = iter_law_convolution(E, ps.m * k)               # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)
    # R2 = build_mod_switching_error_law(ps.q, ps.rq2)    # Rounding2 (in the ciphertext mask part) (ev)
    R2 = build_mod_switching_error_law_AKCN(ps.q, ps.rq2)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)

    # D = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)

    return E, F

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
    dis = ps.q/4-1-ps.q/(2*ps.rq2)
    s = sqrt(s0)
    pr = norm.logsf(dis/s) + 1
    print("err of 3n:")
    print("    = 2^%.2f"% pr)
