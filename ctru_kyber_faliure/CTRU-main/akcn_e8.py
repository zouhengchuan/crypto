from scipy.stats import chi2
from math import sqrt, log, ceil

def ErrorRate(n, l, q, g, sigma, t):
    div_t = 0
    avg_t = 0.
    for i in range(-2**(t-1), 2**(t-1) + 1):
        div_t += (1.*i - 0)**2
    div_t /= 2.**t + 1
    #print t, div_t
    s = sqrt(n*l*sigma**2 * (sigma**2 + div_t) + n*l*sigma**4 + sigma**2)
    dis = (q - 1.) / 2 - sqrt(2.)*(q*1. / g + 1.)
    #print dis, s
    pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(n/16., 2)
    return pr

def ErrorRate2(n, l, q, g, sigma, t):
    div_t = 0
    avg_t = 0.
    avg_t2 = 0.
    for i in range(q):
        i -= (q-1)/2
        temp = i
        temp = temp - q/2**t*round(2**t/q*temp)
        avg_t2 += temp**2/q
        avg_t += temp/q
    div_t = avg_t2 - avg_t
    #print t, div_t
    s = sqrt(n*l*sigma**2 * (sigma**2 + div_t) + n*l*sigma**4 + sigma**2)
    dis = (q - 1.) / 2 - sqrt(2.)*(q*1. / g + 1.) - 1
    #print dis, s
    # pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(n/16., 2)
    pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(n/8., 2)
    return pr

def ErrorRate3(n, l, q, g, sigma, t):
    div_t = 0
    avg_t = 0.
    avg_t2 = 0.
    for i in range(q):
        i -= (q-1)/2
        temp = 2**t/q*i
        temp = temp - round(temp)
        avg_t2 += temp**2/q
        avg_t += temp/q
    div_t = avg_t2 - avg_t
    #print t, div_t
    s = sqrt(1.5*n*l*sigma**2 * (sigma**2 + div_t) + 1.5*n*l*sigma**4 + sigma**2)
    dis = (q - 1.) / 2 - sqrt(2.)*(q*1. / g + 1.) - 1
    #print dis, s
    pr = chi2.logsf((dis/ s)**2, 8.) / log(2) + log(n/8., 2)
    return pr

def build_rounding_law_ciphertext(ps):
    D = {}
    for sigma1 in range(0, ps.q1):    
        temp = 1.*ps.q2/ps.q1*( sigma1 + round(ps.q1/ps.p*k1))
        if temp >= ps.q2:
            temp -= ps.q2
        if temp <= -ps.q2:
            temp += ps.q2                
        epsilon = round( temp ) - temp
        if epsilon >= ps.q2/2:
            epsilon -= ps.q2
        if epsilon <= -ps.q2/2:
            epsilon += ps.q2              
        D[epsilon] = D.get(epsilon,0)+1./ps.q1 * 1./ps.p    
    return D

def bandwidth(n, l, q, g, sigma, t):
    seedLen = 256.
    #print (seedLen + l*n*ceil(log(q/2.**t, 2))) / 8., (l*n*ceil(log(q/2.**t, 2)) + ceil(log(g, 2))*n) / 8.
    return (seedLen + 1.*n*l*ceil(log(q/2.**t, 2)) + 1.*n*l*ceil(log(1.*q, 2)) + ceil((log(g, 2)))*n)/8.

def keySize(n, l, q, g, sigma, t):
    return n*1./2.

n, l, q, g, sigma, t = 1152, 1, 8641, 2**3, sqrt(4.), 10
print (n, "Err rate = ", ErrorRate3(n, l, q, g, sigma, t))
print ("bandwidth = ", bandwidth(n, l, q, g, sigma, t))
print ("keySize = ", keySize(n, l, q, g, sigma, t))