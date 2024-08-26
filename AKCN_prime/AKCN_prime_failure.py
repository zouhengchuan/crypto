from statistics import geometric_mean
from scipy.stats import chi2,  gamma
from math import sqrt, log, ceil, floor

# ErrorRate of prime: x^n - x - 1

# gamma分布逼近
def ErrorRate_1(ps):
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

    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1

    pr = 0
    k_list = []
    theta_list = []
    L1 = []
    L2 = []

    s = sqrt(ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2)
    L1.append(2 * s)
    L2.append(4 * s**2)

    for sig in range(1,ps.n):
        s = (2 * ps.n - sig) * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2
        L1.append(2 * s)
        L2.append(4 * s**2)
    
    for i in range(0,ps.n,8):
        k_list.append(0.5 * sum(L1[i:i+8])**2 / sum(L2[i:i+8]))
        theta_list.append(sum(L2[i:i+8]) / sum(L1[i:i+8]))

    for (k_sum,theta_sum) in zip(k_list,theta_list):
        pr += gamma.sf(dis**2, k_sum, scale=theta_sum)
    pr = log(pr) / log(2)
    print("err of prime (gamma分布逼近):")
    print("    = 2^%.2f"% pr)

# gamma分布逼近,系数按奇偶重排
def ErrorRate_1_prime(ps):
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

    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1

    pr = 0
    k_list = []
    theta_list = []
    L1 = []
    L2 = []

    s = sqrt(ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2)
    L1.append(2 * s)
    L2.append(4 * s**2)

    for sig in range(1,ps.n):
        s = (2 * ps.n - sig) * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2
        L1.append(2 * s)
        L2.append(4 * s**2)
    
    odd_elements = L1[::2]
    even_elements = L1[1::2]
    L1 = odd_elements + even_elements
    odd_elements = L2[::2]
    even_elements = L2[1::2]
    L2 = odd_elements + even_elements

    for i in range(0,ps.n,8):
        k_list.append(0.5 * sum(L1[i:i+8])**2 / sum(L2[i:i+8]))
        theta_list.append(sum(L2[i:i+8]) / sum(L1[i:i+8]))

    for (k_sum,theta_sum) in zip(k_list,theta_list):
        pr += gamma.sf(dis**2, k_sum, scale=theta_sum)
    pr = log(pr) / log(2)
    print("err of prime (gamma分布逼近,系数按奇偶重排):")
    print("    = 2^%.2f"% pr)

# 方差取几何平均
def ErrorRate_2(ps):
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

    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1

    pr = 0
    L = []
    s_list = []
    
    s = sqrt(ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2)
    L.append(s)

    for sig in range(1,ps.n):
        s = sqrt((2 * ps.n - sig ) * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2)
        L.append(s)

    for i in range(0,ps.n,8):
        s_list.append(geometric_mean(L[i:i+8]))
    
    for ss in s_list:
       pr += chi2.sf((dis / ss)**2, 8.)
    pr = log(pr) / log(2)
    print("err of prime (方差几何平均):")
    print("    = 2^%.2f"% pr)

# 方差取几何平均,系数按奇偶重排
def ErrorRate_2_prime(ps):
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

    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1

    pr = 0
    L = []
    s_list = []
    
    s = sqrt(ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2)
    L.append(s)

    for sig in range(1,ps.n):
        s = sqrt((2 * ps.n - sig ) * ps.sigma1**2 * (2 * ps.sigma2**2 + (ps.q/2**ps.t)**2 * div_t) + ps.sigma2**2)
        L.append(s)

    odd_elements = L[::2]
    even_elements = L[1::2]
    L = odd_elements + even_elements

    for i in range(0,ps.n,8):
        s_list.append(geometric_mean(L[i:i+8]))
    
    for ss in s_list:
       pr += chi2.sf((dis / ss)**2, 8.)
    pr = log(pr) / log(2)
    print("err of prime (方差几何平均,系数按奇偶重排):")
    print("    = 2^%.2f"% pr)

#项数取算术平均
def ErrorRate_3(ps):
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

    # standard deviation s, s^2 = 1.5*n*sigma1^2*(2*sigma2^2 + div_t) + sigma2^2
    s = sqrt( (1.5 * ps.n - 1.5)  * ps.sigma1**2 * (2 * ps.sigma2**2 + div_t) + ps.sigma2**2)
    # minimal distance dis = (q - 1)/2 - sqrt(2)*(q/g + 1) - 1
    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1
    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)

        print("err of prime (项数取算术平均):")
        print("    = 2^%.2f"% pr)
    else:
        print("error! dis <= 0")

#项数取最大值
def ErrorRate_4(ps):
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

    # standard deviation s, s^2 = 1.5*n*sigma1^2*(2*sigma2^2 + div_t) + sigma2^2
    s = sqrt((2 * ps.n - 1) * ps.sigma1**2 * (2 * ps.sigma2**2 + div_t) + ps.sigma2**2)
    # minimal distance dis = (q - 1)/2 - sqrt(2)*(q/g + 1) - 1
    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1
    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)

        print("err of prime (项数取最大值):")
        print("    = 2^%.2f"% pr)
    else:
        print("error! dis <= 0")

#项数取最小值
def ErrorRate_5(ps):
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

    # standard deviation s, s^2 = 1.5*n*sigma1^2*(2*sigma2^2 + div_t) + sigma2^2
    s = sqrt( ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + div_t) + ps.sigma2**2)
    # minimal distance dis = (q - 1)/2 - sqrt(2)*(q/g + 1) - 1
    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1
    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)

        print("err of prime (项数取最大值):")
        print("    = 2^%.2f"% pr)
    else:
        print("error! dis <= 0")


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