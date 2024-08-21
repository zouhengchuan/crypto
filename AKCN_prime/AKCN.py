from math import log
import operator as op
from math import factorial as fac
from math import sqrt, log, ceil
from proba_util import *
from AKCN_failure import *
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet


class AKCN_ParameterSet:
    def __init__(self, n, q, g, t, eta_1, eta_2=None):
        if eta_2 is None:
            eta_2 = eta_1
        self.n = n
        self.q = q
        self.g = g
        self.t = t
        self.eta1 = eta_1
        self.eta2 = eta_2
        # standard deviation used in LWE attack
        self.sigma = sqrt(sqrt(eta_1*eta_2)/2)
        # standard deviation of CBD(eta_1) and CBD(eta_2)
        self.sigma1 = sqrt(var_of_law(build_centered_binomial_law(eta_1)))
        self.sigma2 = sqrt(var_of_law(build_centered_binomial_law(eta_2)))
        

def AKCN_to_MLWE(ps):
    return MLWEParameterSet(ps.n, ps.q, ps.sigma)

def summarize(ps):
    MLWE_summarize_attacks(AKCN_to_MLWE(ps))
    ErrorRate(ps)
    ErrorRate_new(ps)
    f = ps_probability(ps)
    print ("failure origin: %.1f = 2^%.5f"%(f, log(f + 2.**(-300))/log(2)))
    

if __name__ == "__main__":

    # ParameterSet(n, q, g, t, eta_1, eta_2)
    L = []
    # for k in range(3,6):
    #     for t in range(8,12):
    #         ps_test = AKCN_ParameterSet(1152, 3457, 2**k, t, 2, 2)
    #         L.append(ps_test)
    ps_test = AKCN_ParameterSet(768, 3457, 2**3, 10, 2, 2)
    L.append(ps_test)

    for ps_temp in L:
        print("**************************************")
        print("AKCN: n = %d, q = %d, g = %d, t = %d, eta_1 = %d, eta_2 = %d\n"%(ps_temp.n, ps_temp.q, ps_temp.g, ps_temp.t, ps_temp.eta1, ps_temp.eta2))  
        bandwidth(ps_temp)
        summarize(ps_temp)
        print()