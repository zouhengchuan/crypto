from scipy.stats import chi2
from math import sqrt, log, ceil, floor
from proba_util import*

# ErrorRate of x^n - x^{n/2} + 1
def ErrorRate(ps):
    # variance of (q/2^t)*epsilon = (2^t/q)*u - round(2^t/q*u)
    # u is uniform in (-(q-1)/2,(q-1)/2), round(x) = floor(x + 1/2)
    div_t = 0
    avg_t = 0
    for i in range(-((ps.q-1)//2),(ps.q-1)//2+1):
        temp = 2**ps.t / ps.q * i - floor(2**ps.t / ps.q * i + 1/2)
        avg_t += temp / ps.q
    for i in range(-((ps.q-1)//2),(ps.q-1)//2+1):
        temp = 2**ps.t / ps.q * i - floor(2**ps.t / ps.q * i + 1/2)
        div_t += (temp-avg_t)**2 / ps.q

    # standard deviation s, s^2 = 1.5*n*sigma1^2*(2*sigma2^2 + div_t) + sigma2^2
    s = sqrt(1.5 * ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2)
    # minimal distance dis = (q - 1)/2 - sqrt(2)*(q/g + 1) - 1
    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1
    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)

        print("err of 3n:")
        print("    = 2^%.2f"% pr," (dis/s)^2 = %.2f"%(dis/ s)**2," div_t = %.2f"%div_t," dis = %.2f"%dis," s = %.2f"%s)
    else:
        print("error! dis <= 0")

def ErrorRate_new(ps):

    div_t = var_of_law(build_mod_switching_error_law(ps.q, 2**ps.t))

    # standard deviation s, s^2 = 1.5*n*sigma1^2*(2*sigma2^2 + div_t) + sigma2^2
    s = sqrt(1.5 * ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + div_t) + ps.sigma2**2)
    # minimal distance dis = (q - 1)/2 - sqrt(2)*(q/g + 1) - 1
    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.)
    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)

        print("err of 3n:")
        print("    = 2^%.2f"% pr," (dis/s)^2 = %.2f"%(dis/ s)**2," div_t = %.2f"%div_t," dis = %.2f"%dis," s = %.2f"%s)
    else:
        print("error! dis <= 0")

# ErrorRate of x^n + 1
def ErrorRate_2(ps):
    # variance of (q/2^t)*epsilon = u - (q/2^t)*round(2^t/q*u)
    # u is uniform in (-(q-1)/2,(q-1)/2), round(x) = floor(x + 1/2)
    div_t = 0
    avg_t = 0
    for i in range(-((ps.q-1)//2),(ps.q-1)//2+1):
        temp = i - ps.q / 2**ps.t * floor(2**ps.t / ps.q * i + 1/2)
        avg_t += temp / ps.q
    for i in range(-((ps.q-1)//2),(ps.q-1)//2+1):
        temp = i - ps.q / 2**ps.t * floor(2**ps.t / ps.q * i + 1/2)
        div_t += (temp-avg_t)**2 / ps.q

    # standard deviation s, s^2 = n*sigma1^2*(2*sigma2^2 + div_t) + sigma2^2
    s = sqrt(ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + div_t) + ps.sigma2**2)
    # minimal distance dis = (q - 1)/2 - sqrt(2)*(q/g + 1) - 1
    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1
    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)
        print("err of 2n:")
        print("    = 2^%.2f"% pr)
    else:
        print("error! dis <= 0")

# distribution
def ps_distribution(ps):
    """ construct the final error distribution in our encryption scheme
    :param ps: parameter set (ParameterSet)
    """
    chis = build_centered_binomial_law(ps.eta1)        # LWE error law for the key (x)
    chie = build_centered_binomial_law(ps.eta2)        # (e)

    Rc = build_mod_switching_error_law(ps.q, 2**ps.t)   # rounding error first ciphertext (eu)
    chiRe = law_convolution(chie, Rc)                   # LWE + rounding error ciphertext  (e1 + eu)

    B1 = law_product(chis, chie)
    B2 = law_product(chis, chiRe)                       # (s*e1 + s*eu)

    C1 = iter_law_convolution(B1, ps.n * 3 // 2)          # m*n*(e*r + et*r) prime*2
    C2 = iter_law_convolution(B2, ps.n * 3 // 2)          # m*n*(s*e1 + s*eu) 

    C=law_convolution(C1, C2)                           # m*n*(e*r+et*r) + m*n*(s*e1+s*eu)

    R2 = build_mod_switching_error_law(ps.q, ps.g)    # Rounding2 (in the ciphertext mask part) (ev)
    F = law_convolution(R2, chie)                       # LWE+Rounding2 error (ev + e2)
    D1 = law_convolution(C, F)                           # Final error (m*n*(e*r+et*r) + m*n*(s*e1+s*eu) + ev + e2)

    D2 = build_law_square(D1, ps.q)
    # D = iter_law_convolution(D2,8)
    D = law_convolution(D2,D2)
    D = law_convolution(D,D)
    D = law_convolution(D,D)
    return D

def ps_probability(ps):
    F = ps_distribution(ps)
    #proba = tail_probability(F, ps.q/4)
    proba = tail_probability_frac(F, 0, (ps.q-1)**2/4)
    return ps.n*proba/8

def bandwidth(ps):
    seedLen = 256
    pk = seedLen//8 + ceil(1.*ps.n*ceil(log(ps.q,2))/8)
    ct = ceil(ps.n*ps.t/8) + ceil(ps.n*ceil(log(ps.g,2))/8)
    print('|pk| = %d, |ct| = %d, bandwidth = %d\n'%(pk, ct, (pk+ct)))
def keySize(ps):
    return ps.n*1./2.