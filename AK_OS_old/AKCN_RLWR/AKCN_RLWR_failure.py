from statistics import geometric_mean
from scipy.stats import chi2,  gamma
from math import sqrt, log, ceil, floor
from proba_util import*
import matplotlib.pyplot as plt

def distribution(ps):
    A = build_rq_law(2**ps.eq)
    X1 = build_centered_binomial_law(ps.eta)
    # X2 = build_centered_binomial_law(ps.eta)

    AX1 = law_product_over_non_power_of_2(A, X1, 2**ps.eq)
    AX1_3n = iter_law_convolution_q(AX1, ps.n//2, 2**ps.eq)
    AX1_p = err_law(AX1_3n, 2**ps.eq, 2**ps.ep)

    # AX2 = law_product_over_non_power_of_2(A, X2, 2**ps.eq)
    # AX2_3n = iter_law_convolution_q(AX2, ps.n//2, 2**ps.eq)
    # AX2_p = err_law(AX2_3n, 2**ps.eq, 2**ps.ep)

    # err1 = law_product_over_non_power_of_2(AX2_p, X1, 2**ps.eq)
    # err1_3n = iter_law_convolution_q(err1, ps.n//2, 2**ps.eq)
    # err2 = law_product_over_non_power_of_2(AX1_p, X2, 2**ps.eq)
    # err2_3n = iter_law_convolution_q(err2, ps.n//2, 2**ps.eq)

    err = law_product_over_non_power_of_2(AX1_p, X1, 2**ps.eq)
    err_3n = iter_law_convolution_q(err, ps.n//2, 2**ps.eq)

    sig = law_nconvolution_q(err_3n, err_3n, 2**ps.eq)
    sig_sq = build_law_square(sig, 2**ps.eq)

    dis = (2**ps.ep / 2 - sqrt(2) * 2**(ps.ep - ps.eg))**2
    sig_8 = iter_law_convolution_re(sig_sq, 8, dis)

    # fig, axs = plt.subplots(2, 1, figsize=(12, 7))
    # axs[0, 0].plot(err_3n.keys(), err_3n.values(), marker='o', markersize=1, linestyle='', color = "red")
    # axs[0, 0].set_title("err_3n")

    # axs[1, 0].plot(sig_8.keys(), sig_8.values(), marker='o', markersize=1, linestyle='', color = "blue")
    # axs[1, 0].set_title("sig_8")

    # plt.suptitle('(q, p) = (%d, %d)'%(ps.ep, ps.eq), fontsize=16)
    # plt.tight_layout()
    # plt.show()

    return sig_8

def prob(ps):
    D = distribution(ps)
    dis = (2**ps.ep / 2 - sqrt(2) * 2**(ps.ep - ps.eg))**2
    proba = tail_probability_frac(D, 0, dis)
    pp =  ps.n * proba / 8
    logpp = log(pp, 2)
    print("err:")
    print("    = 2^%.2f"% logpp)

# ErrorRate of AKCN-RLWR: x^n - x^(n/2) + 1


#项数取最大值
def ErrorRate(ps):
    # standard deviation s, s^2 = 3 * n * (p/q)^2 * sigma1^2 * sigma2^2
    s = sqrt(3 * ps.n * 2**(2 * (ps.ep - ps.eq)) * ps.sigma1**2 * ps.sigma2**2)
    # minimal distance dis = p/2 - sqrt(2) * p/g 
    dis = (2**ps.ep) / 2 - sqrt(2.) * 2**(ps.ep - ps.eg) 

    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)

        print("err of AKCN-RLWR (项数取最大值):")
        print("    = 2^%.2f"% pr)
    else:
        print("error! dis <= 0")


def ErrorRate_test(ps):
    # minimal distance dis = p/2 - sqrt(2) * p/g 
    dis = (2**ps.ep) / 2 - sqrt(2.) * 2**(ps.ep - ps.eg) 
    
    pr = 0
    k_list = []
    theta_list = []
    L1 = []
    L2 = []

    for sig in range(ps.n//2):
        s = (1.5 * ps.n - sig - 1) * 2 * 2**(2 * (ps.ep - ps.eq)) * ps.sigma1**2 * ps.sigma2**2
        L1.append(2 * s)
        L2.append(4 * s**2)
    for sig in range(ps.n//2):
        s = (1.5 * ps.n) * 2 * 2**(2 * (ps.ep - ps.eq)) * ps.sigma1**2 * ps.sigma2**2
        L1.append(2 * s)
        L2.append(4 * s**2)
    for i in range(0,ps.n,8):
        k_list.append(0.5 * sum(L1[i:i+8])**2 / sum(L2[i:i+8]))
        theta_list.append(sum(L2[i:i+8]) / sum(L1[i:i+8]))

    for (k_sum,theta_sum) in zip(k_list,theta_list):
        pr += gamma.sf(dis**2, k_sum, scale=theta_sum)
    
    pr = log(pr) / log(2)
    print("err of test(gamma分布逼近):")
    print("    = 2^%.2f"% pr)



def bandwidth(ps):
    seedLen = ps.n/16.
    pk = seedLen + ceil( ps.n * ps.ep/8)
    ct = ceil(ps.n*ps.ep/8) + ceil(ps.n*ps.eg/8)
    print('|pk| = %d, |ct| = %d, bandwidth = %d\n'%(pk, ct, (pk+ct)))

def bandwidth_K320(ps):
    seedLen = 40.
    pk = seedLen + ceil(1.*ps.n*ceil(log(ps.q,2))/8)
    ct = ceil(ps.n*ps.t/8) + ceil(ps.n*ceil(log(ps.g,2))/8)
    print('|pk| = %d, |ct| = %d, bandwidth = %d\n'%(pk, ct, (pk+ct)))

def bandwidth_K512(ps):
    seedLen = 64.
    pk = seedLen + ceil(1.*ps.n*ceil(log(ps.q,2))/8)
    ct = ceil(ps.n*ps.t/8) + ceil(ps.n*ceil(log(ps.g,2))/8)
    print('|pk| = %d, |ct| = %d, bandwidth = %d\n'%(pk, ct, (pk+ct)))

def keySize(ps):
    return ps.n*1./2.