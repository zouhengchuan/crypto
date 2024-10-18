from math import factorial as fac
from math import log, ceil, erf, sqrt, floor


''' Centered binomial distribution, parameter eta = 1 '''
def build_cbd1_law():
    D = build_centered_binomial_law(1)
    return D

def build_cbd2_law():
    D = build_centered_binomial_law(2)
    return D

def build_cbd3_law():
    D = build_centered_binomial_law(3)
    return D

''' Centered binomial distribution, parameter eta = 4 '''
def build_cbd4_law():
    D = build_centered_binomial_law(4)
    return D

def build_cbd5_law():
    D = build_centered_binomial_law(5)
    return D

''' Calculate the variance of distribution D '''
def var_of_law(D):
    div_t = 0
    avg_t = 0.
    for d in D:
        avg_t += d
    avg_t /= len(D)
    for d in D:
        div_t += D[d]*(d - avg_t)**2
    return div_t
 

def gaussian_center_weight(sigma, t):
    """ Weight of the gaussian of std deviation s, on the interval [-t, t]
    :param x: (float)
    :param y: (float)
    :returns: erf( t / (sigma*\sqrt 2) )
    """
    return erf(t / (sigma * sqrt(2.)))


def binomial(x, y):
    """ Binomial coefficient
    :param x: (integer)
    :param y: (integer)
    :returns: y choose x
    """
    try:
        binom = fac(x) // fac(y) // fac(x - y)
    except ValueError:
        binom = 0
    return binom


def centered_binomial_pdf(k, x):
    """ Probability density function of the centered binomial law of param k at x
    :param k: (integer)
    :param x: (integer)
    :returns: p_k(x)
    """
    return binomial(2*k, x+k) / 2.**(2*k)


def build_centered_binomial_law(k):
    """ Construct the binomial law as a dictionnary
    :param k: (integer)
    :param x: (integer)
    :returns: A dictionnary {x:p_k(x) for x in {-k..k}}
    """
    D = {}
    for i in range(-k, k+1):
        D[i] = centered_binomial_pdf(k, i)
    return D

def build_rounding_law_rlwr(eq,ep):
    D = {}
    t = eq - ep - 1
    for u in range(-2**t, 2**t):    
        epsilon = u
        D[epsilon] = D.get(epsilon,0)+2**(ep - eq)
    return D 

def build_law_square(A, q):
    C = {}
    for a in A:
        c = mod_centered(a, q)**2
        C[c] = C.get(c, 0) + A[a]
    return C

def build_rq_law(q):
    D = {}
    for i in range(q):
        D[i] = 1./q
    return D

def err_law(AX, q, p):
    D = {}
    for i in AX:
        err = p/q * i - floor(p/q * i + 1./2)
        D[err] = D.get(err,0) + AX[i]
    return D

def mod_switch(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(round(1.* rq * x / q) % rq)


def mod_centered(x, q):
    """ reduction mod q, centered (ie represented in -q/2 .. q/2)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    """
    a = x % q
    if a < q/2:
        return a
    return a - q


def build_mod_switching_error_law(q, rq):
    """ Construct Error law: law of the difference introduced by switching from and back a uniform value mod q
    :param q: original modulus (integer)
    :param rq: intermediate modulus (integer)
    """
    D = {}
    V = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch(y, rq, q)
        d = mod_centered(x - z, q)
        D[d] = D.get(d, 0) + 1./q
        V[y] = V.get(y, 0) + 1

    return D


def law_convolution(A, B):
    """ Construct the convolution of two laws (sum of independent variables from two input laws)
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """

    C = {}
    for a in A:
        for b in B:
            c = a+b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def law_convolution_re(A, B, dis):
    C = {}
    for a in A:
        for b in B:
            c = floor((a+b)/5) * 5
            if (c > dis):
                c = dis + 1
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def law_convolution_q(A, B, q):
    C = {}
    for a in A:
        for b in B:
            c = (a+b) % q
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def law_nconvolution_q(A, B, q):
    C = {}
    for a in A:
        for b in B:
            c = (a-b) % q
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C

def law_product(A, B):
    """ Construct the law of the product of independent variables from two input laws
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """
    C = {}
    for a in A:
        for b in B:
            c = a*b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C
 
def law_product_over_non_power_of_2(A, B, q):
    """ Construct the law of the non-power-of-2 product of independent variables from two input laws
    :param A: first input law (dictionnary)
    :param B: second input law (dictionnary)
    """
    C = {}
    for a in A:
        for a_prime in A:
            for b in B:
                for b_prime in B:
                    c = (a * b + b_prime * (a + a_prime)) % q
                    C[c] = C.get(c, 0) + A[a] * B[b] * A[a_prime] * B[b_prime]
    return C


def clean_dist(A):
    """ Clean a distribution to accelerate further computation (drop element of the support with proba less than 2^-300)
    :param A: input law (dictionnary)
    """
    B = {}
    for (x, y) in A.items():
        if y>2**(-300):
            B[x] = y
    return B


def iter_law_convolution(A, i):
    """ compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution(D, D)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution(D, A)
            D = clean_dist(D)
    return D

def iter_law_convolution_re(A, i, dis):
    """ compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution_re(D, D, dis)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution_re(D, A, dis)
            D = clean_dist(D)
    return D

def iter_law_convolution_q(A, i, q):
    """ compute the -ith forld convolution of a distribution (using double-and-add)
    :param A: first input law (dictionnary)
    :param i: (integer)
    """
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution_q(D, D, q)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution_q(D, A, q)
            D = clean_dist(D)
    return D

def tail_probability(D, t):
    '''
    Probability that an drawn from D is strictly greater than t in absolute value
    :param D: Law (Dictionnary)
    :param t: tail parameter (integer)
    '''
    s = 0
    ma = max(D.keys())
    if t >= ma:
        return 0
    for i in reversed(range(int(ceil(t)), ma)):  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        s += D.get(i, 0) + D.get(-i, 0)
    return s


def tail_probability_frac(D, begin, end):
    '''
    Probability that an drawn from D is strictly greater than t in absolute value
    :param D: Law (Dictionnary)
    :param t: tail parameter (integer)
    '''
    s = 0
    ma = max(D.keys())
    mi = min(D.keys())
    if (begin < mi) and (end > ma):
        return 0
    for d in D:  # Summing in reverse for better numerical precision (assuming tails are decreasing)
        if (d < begin) or (d >= end):
            s += D.get(d, 0)
    return s