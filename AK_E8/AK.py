from math import log
import operator as op
from math import factorial as fac
from math import sqrt, log, ceil
from proba_util import *
from AK_failure import *
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet


class AKCN_ParameterSet:
    def __init__(self, n, q, g, t_1, t_2, eta_1, eta_2=None):
        if eta_2 is None:
            eta_2 = eta_1
        self.n = n
        self.q = q
        self.g = g
        self.t1 = t_1
        self.t2 = t_2
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
    

if __name__ == "__main__":

    # ParameterSet(n, q, g, t_1, t_2, eta_1, eta_2)
    L = []
    for k in range(3,5):
        for t_1 in range(1,4):
            for t_2 in range (t_1,t_1+2):
                ps_test = AKCN_ParameterSet(768, 3457, 2**k, t_1, t_2, 2, 2)
                L.append(ps_test)

    # L = [AKCN_ParameterSet(768, 3457, 2**4, 2, 2, 3, 3)]
    for ps_temp in L:
        print("**************************************")
        print("AKCN: n = %d, q = %d, g = %d, t_1 = %d, t_2 = %d, eta_1 = %d, eta_2 = %d\n"%(ps_temp.n, ps_temp.q, ps_temp.g, ps_temp.t1, ps_temp.t2, ps_temp.eta1, ps_temp.eta2))  
        bandwidth(ps_temp)
        summarize(ps_temp)
        print()