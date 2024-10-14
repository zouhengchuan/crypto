from math import log
import operator as op
from math import factorial as fac
from math import sqrt, log, ceil
from proba_util import *
from AKCN_RLWR_failure import *
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet


class AKCN_ParameterSet:
    def __init__(self, n, eq, ep, eg, eta):
        self.n = n
        self.eq = eq
        self.ep = ep
        self.eg = eg
        self.eta = eta
        self.sigma1 = sqrt(var_of_law(build_centered_binomial_law(eta)))
        self.sigma2 = sqrt(var_of_law(build_rounding_law_rlwr(eq,ep)))

        # standard deviation used in LWE attack
        self.sigma = sqrt(self.sigma1 * self.sigma2)
        
        

'''
def AKCN_to_MLWE(ps):

    # Check whether ciphertext error variance after rounding is larger than secret key error variance
    Rc = build_mod_switching_error_law(ps.q, 2**ps.t)
    var_rounding = sum([i*i*Rc[i] for i in Rc.keys()])

    if ps.eta2/2. + var_rounding < ps.eta1/2.:
        raise "The security of the ciphertext MLWE may not be stronger than the one of the public key MLWE"    

    return MLWEParameterSet(ps.n, ps.q, ps.sigma) 
'''    

def AKCN_to_MLWE(ps):
    return MLWEParameterSet(ps.n, 2**ps.eq, ps.sigma)

def summarize(ps):
    MLWE_summarize_attacks(AKCN_to_MLWE(ps))
    # ErrorRate(ps)
    ErrorRate_test(ps)
    

if __name__ == "__main__":

    # ParameterSet(n, eq, ep, eg, eta)
    L = []
    for k in range(10,11):
       for t in range(2,10):
           ps_test = AKCN_ParameterSet(768, 13, 10, t, 4)
           L.append(ps_test)
    
    # ps_test1 = AKCN_ParameterSet(768, 11, 9, 4, 2)
    # ps_test1 = AKCN_ParameterSet(768, 11, 9, 4, 2)
    # L.append(ps_test1)

    # ps_test2 = AKCN_ParameterSet(1152, 13, 10, 3, 2)
    # L.append(ps_test2)
 
    for ps_temp in L:
        print("**************************************")
        print("AKCN_RLWR: n = %d, eq = %d, ep = %d, eg = %d, eta = %d\n"%(ps_temp.n, ps_temp.eq, ps_temp.ep, ps_temp.eg, ps_temp.eta))  
        bandwidth(ps_temp)
        summarize(ps_temp)
        print("sigma1&2=",ps_temp.sigma1,ps_temp.sigma2)
        print()
    
    # L_K320 = []
    #for k in range(3,4):
    #    for t in range(8,9):
    #        ps_test = AKCN_ParameterSet(761, 4091, 2**k, t, 2, 3)
    #       L_K320.append(ps_test)

    #ps_test = AKCN_ParameterSet(769, 6599, 2**13, 9, 3, 3)
    #L_K320.append(ps_test)

    # ps_test1_K320 = AKCN_ParameterSet(751, 3823, 2**4, 10, 3, 3)
    # L_K320.append(ps_test1_K320)

    #ps_test2_K320 = AKCN_ParameterSet(757, 3727, 2**4, 10, 2, 3)
    #L_K320.append(ps_test2_K320)

    # for ps_temp in L_K320:
    #     print("**************************************")
    #     print("AKCN_prime: n = %d, q = %d, g = %d, t = %d, eta_1 = %d, eta_2 = %d\n"%(ps_temp.n, ps_temp.q, ps_temp.g, ps_temp.t, ps_temp.eta1, ps_temp.eta2))  
    #     bandwidth(ps_temp)
    #     summarize(ps_temp)
    #     print()

    # ParameterSet(n, q, g, t, eta_1, eta_2)
    #L_K512 = []
    #for k in range(3,4):
    #    for t in range(8,9):
    #        ps_test = AKCN_ParameterSet(761, 4091, 2**k, t, 2, 3)
    #       L_K512.append(ps_test)

    #ps_test = AKCN_ParameterSet(769, 6599, 2**13, 9, 3, 3)
    #L_K512.append(ps_test)

    #for ps_temp in L_K512:
        #print("**************************************")
        #print("AKCN_prime: n = %d, q = %d, g = %d, t = %d, eta_1 = %d, eta_2 = %d\n"%(ps_temp.n, ps_temp.q, ps_temp.g, ps_temp.t, ps_temp.eta1, ps_temp.eta2))  
        #bandwidth_K512(ps_temp)
        #summarize(ps_temp)
        #print()