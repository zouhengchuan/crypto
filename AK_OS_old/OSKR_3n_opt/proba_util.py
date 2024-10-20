from math import factorial as fac
from math import log, ceil, floor, erf, sqrt

def center_floor(x):
    """
    the x in operation x%y will be coverted into integer
    """
    y = x
    if y <= 0:
        y = ceil(y)
    if y > 0:
        y = floor(y)
    return y

def mod_q(x,q):
    """
    the x in operation x%y will be coverted into integer
    """
    y = x
    while (y <= 0):
        y = y + q
    while (y >= q):
        y = y - q
    return y


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


def mod_switch(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(floor(1.* rq * x / q + 0.5) % rq)


def mod_switch_floor(x, q, rq):
    """ Modulus switching (rounding to a different discretization of the Torus)
    :param x: value to round (integer)
    :param q: input modulus (integer)
    :param rq: output modulus (integer)
    """
    return int(center_floor(1.* rq * x / q) % rq)


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
    """ Construct Error law: law of the difference introduced by con and rec a uniform value mod q
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

def build_mod_switching_error_law_floor(q, rq):
    """ Construct Error law: law of the difference introduced by con and rec a uniform value mod q
    :param q: original modulus (integer)
    :param rq: intermediate modulus (integer)
    """
    D = {}
    V = {}
    for x in range(q):
        y = mod_switch(x, q, rq)
        z = mod_switch_floor(y, rq, q)
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


def clean_dist(A):
    """ Clean a distribution to accelerate further computation (drop element of the support with proba less than 2^-300)
    :param A: input law (dictionnary)
    """
    B = {}
    for (x, y) in A.items():
        if y>2**(-300):
            B[x] = y
    return B


def iter_law_product_3n(A, X, n, q, k):
    AX_1 = law_product_3n(A, X)
    AX_2 = law_product_3n_2(A, X)
    AX_11 = iter_law_convolution_q(AX_1, (n//2 - k), q)
    AX_22 = iter_law_convolution_q(AX_2, k, q)
    AX = law_convolution_q(AX_11, AX_22, q)
    return AX


def law_convolution_q(A, B, q):
    C = {}
    for a in A:
        for b in B:
            # c = mod_centered(a + b, q)
            c = a + b
            C[c] = C.get(c, 0) + A[a] * B[b]
    return C


def iter_law_convolution_q(A, i, q):
    D = {0: 1.0}
    i_bin = bin(i)[2:]  # binary representation of n
    for ch in i_bin:
        D = law_convolution_q(D, D, q)
        D = clean_dist(D)
        if ch == '1':
            D = law_convolution_q(D, A, q)
            D = clean_dist(D)
    return D


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
        if (d < begin) or (d > end):
            s += D.get(d, 0)
    return s
