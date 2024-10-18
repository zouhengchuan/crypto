from scipy.stats import gamma
from math import sqrt, log, ceil, floor
from proba_util import*

# ErrorRate of AKCN-E8-Hybrid-prime: x^n - x - 1

def ErrorRate_Hybrid_prime(ps):

    # minimal distance dis = p/2 - sqrt(2) * p/g 

    dis = ps.p / 2 - sqrt(2.) * (ps.p / ps.g) 
    
    pr = 0
    k_list = []
    theta_list = []
    L1 = []
    L2 = []

    # {ar}_p = q/p * round(p/q * a * r ) - a * r
    # p/q * e * r - p/q * {ar}_p * s
    
    s = ps.n * (ps.p / ps.q)**2 * ( ps.sigmaE2**2 * ps.sigmaR1**2 + ps.sigmaE1**2 * ps.sigmaR2**2 )
    L1.append(2 * s)
    L2.append(4 * s**2)

    for sig in range(1,ps.n):
        s = (2 * ps.n - sig) * (ps.p / ps.q)**2 * ( ps.sigmaE2**2 * ps.sigmaR1**2 + ps.sigmaE1**2 * ps.sigmaR2**2 )
        L1.append(2 * s)
        L2.append(4 * s**2)

    odd_elements = L1[::2]
    even_elements = L1[1::2]
    L1 = odd_elements + even_elements
    odd_elements = L2[::2]
    even_elements = L2[1::2]
    L2 = odd_elements + even_elements

    # gamma

    ue = 16 * floor(ps.n / 16)
    for i in range(0,ue,8):
        k_list.append(0.5 * sum(L1[i:i+8])**2 / sum(L2[i:i+8]))
        theta_list.append(sum(L2[i:i+8]) / sum(L1[i:i+8]))

    for (k_sum,theta_sum) in zip(k_list,theta_list):
        pr += gamma.sf(dis**2, k_sum, scale=theta_sum)
    
    pr = log(pr) / log(2)
    print("err of AKCN-E8-Hybrid-prime:")
    print("    = 2^%.2f"% pr)


def bandwidth(ps):
    seedLen = floor(ps.n/16)
    pk = seedLen + ceil( ps.n * log(ps.q,2)/8)
    ct = ceil(ps.n*log(ps.p,2)/8) + ceil(16 * seedLen * log(ps.g,2)/8)
    print('|K| = %d, |pk| = %d, |ct| = %d, bandwidth = %d\n'%(seedLen, pk, ct, (pk+ct)))

def bandwidth_K320(ps):
    seedLen = 40.
    pk = seedLen + ceil( ps.n*ceil(log(ps.q,2))/8 )
    ct = ceil(ps.n*log(ps.p,2)/8) + ceil(640*log(ps.g,2)/8)
    print('|K| = %d, |pk| = %d, |ct| = %d, bandwidth = %d\n'%(seedLen, pk, ct, (pk+ct)))

def bandwidth_K512(ps):
    seedLen = 64.
    pk = seedLen + ceil( ps.n*ceil(log(ps.q,2))/8 )
    ct = ceil(ps.n*log(ps.p,2)/8) + ceil(1024*log(ps.g,2)/8)
    print('|K| = %d, |pk| = %d, |ct| = %d, bandwidth = %d\n'%(seedLen, pk, ct, (pk+ct)))

def keySize(ps):
    return ps.n*1./2.