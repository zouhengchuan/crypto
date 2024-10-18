from scipy.stats import gamma
from math import sqrt, log, ceil
from proba_util import*

# ErrorRate of AKCN-E8-Hybrid: x^n - x^(n/2) + 1

def ErrorRate_Hybrid(ps):

    # minimal distance dis = p/2 - sqrt(2) * p/g 

    dis = ps.p / 2 - sqrt(2.) * (ps.p / ps.g) 
    
    pr = 0
    k_list = []
    theta_list = []
    L1 = []
    L2 = []
    
    # {ar}_p = q/p * round(p/q * a * r ) - a * r
    # p/q * e * r - p/q * {ar}_p * s

    for sig in range(ps.n//2):
        s = (1.5 * ps.n - sig - 1) * (ps.p / ps.q)**2 * ( ps.sigmaE2**2 * ps.sigmaR1**2 + ps.sigmaE1**2 * ps.sigmaR2**2 )
        L1.append(2 * s)
        L2.append(4 * s**2)
    for sig in range(ps.n//2):
        s = (1.5 * ps.n) * (ps.p / ps.q)**2 * ( ps.sigmaE2**2 * ps.sigmaR1**2 + ps.sigmaE1**2 * ps.sigmaR2**2 )
        L1.append(2 * s)
        L2.append(4 * s**2)

    # gamma

    for i in range(0,ps.n,8):
        k_list.append(0.5 * sum(L1[i:i+8])**2 / sum(L2[i:i+8]))
        theta_list.append(sum(L2[i:i+8]) / sum(L1[i:i+8]))

    for (k_sum,theta_sum) in zip(k_list,theta_list):
        pr += gamma.sf(dis**2, k_sum, scale=theta_sum)
    
    pr = log(pr) / log(2)
    print("err of AKCN-E8-Hybrid:")
    print("    = 2^%.2f"% pr)


def bandwidth(ps):
    seedLen = ps.n/16.
    pk = seedLen + ceil( ps.n * log(ps.q,2)/8)
    ct = ceil(ps.n*log(ps.p,2)/8) + ceil(ps.n*log(ps.g,2)/8)
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