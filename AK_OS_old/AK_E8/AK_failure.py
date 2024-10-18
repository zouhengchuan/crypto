from scipy.stats import chi2
from math import sqrt, log, ceil, floor

def var_epsilon(q,t):
    # variance of 2^t*epsilon_t = 2^t*round(u/2^t) - u
    # u is uniform in (-(q-1)/2,(q-1)/2), round(x) = floor(x + 1/2)
    div_t = 0
    avg_t = 0
    for u in range(-((q-1)//2),(q-1)//2+1):
        temp = 2**t * floor(u / 2**t + 1/2) - u
        avg_t += temp / q
    for u in range(-((q-1)//2),(q-1)//2+1):
        temp = 2**t * floor(u / 2**t + 1/2) - u
        div_t += (temp-avg_t)**2 / q
    return div_t

# ErrorRate of x^n - x^{n/2} + 1
def ErrorRate(ps):
    div_t1 = var_epsilon(ps.q,ps.t1)
    div_t2 = var_epsilon(ps.q,ps.t2)

    # standard deviation s, s^2 = 1.5*n*sigma1^2*(2*sigma2^2 + div_t1 + div_t2) + sigma2^2
    s = sqrt(1.5 * ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + div_t1 + div_t2) + ps.sigma2**2)

    # minimal distance dis = (q - 1)/2 - sqrt(2)*(q/g + 1) - 1
    dis = (ps.q - 1.) / 2 - sqrt(2.)*(ps.q*1. / ps.g + 1.) - 1

    # error rate in each block is delta_0 = Pr(d from chi(8) | sqrt(d) > dis/s)
    # then error rate is delta = 1 - (1 - delta_0)^(n/8), approximate to n/8*delta_0
    # so log(delta,2) = logsf(delta0)/log(2) + log(n/8,2)
    if dis > 0:
        pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(ps.n/8., 2)

        print("err of 3n:")
        print("    = 2^%.2f"% pr)
    else:
        print("error! dis <= 0")

# ErrorRate of x^n + 1
def ErrorRate_2(ps):
    div_t1 = var_epsilon(ps.q,ps.t1)
    div_t2 = var_epsilon(ps.q,ps.t2)

    # standard deviation s, s^2 = n*sigma1^2*(2*sigma2^2 + div_t1 + div_t2) + sigma2^2
    s = sqrt(ps.n * ps.sigma1**2 * (2 * ps.sigma2**2 + div_t1 + div_t2) + ps.sigma2**2)

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

def bandwidth(ps):
    seedLen = 256.
    pk = seedLen + ceil(1.*ps.n*ceil(log(ps.q,2)-ps.t1)/8)
    ct = ceil(ps.n*ceil(log(ps.q,2)-ps.t2)/8) + ceil(ps.n*ceil(log(ps.g,2))/8)
    print('|pk| = %d, |ct| = %d, bandwidth = %d\n'%(pk, ct, (pk+ct)))
def keySize(ps):
    return ps.n*1./2.