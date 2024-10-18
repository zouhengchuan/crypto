from math import sqrt, log
from proba_util import *
from AKCN_RLWR_failure import *
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet


class AKCN_ParameterSet:
    def __init__(self, n, q, p, g, eta):
        self.n = n
        self.q = q
        self.p = p
        self.g = g
        self.eta = eta
        self.sigma1 = sqrt(var_of_law(build_centered_binomial_law(eta)))
        self.sigma2 = sqrt(var_of_law(build_rounding_law_rlwr_powerof2(q,p)))

        # standard deviation used in LWE attack
        self.sigma = sqrt(self.sigma1 * self.sigma2)
        
        

def AKCN_to_MLWE(ps):
    return MLWEParameterSet(ps.n, ps.q, ps.sigma)

def summarize(ps):
    MLWE_summarize_attacks(AKCN_to_MLWE(ps))
    ErrorRate_RLWR(ps)
    

if __name__ == "__main__":

    # ParameterSet(n, q, p, g, eta)
    # |K| = n/16
    L = []

    #for k in range(,):
    #    for t in range(,):
    #        ps_test = AKCN_ParameterSet(, , , , , )
    #        L.append(ps_test)

    ps_test1_1 = AKCN_ParameterSet(768, 2**12, 2**9, 2**5, 2)
    L.append(ps_test1_1)
    ps_test1_2 = AKCN_ParameterSet(768, 2**11, 2**9, 2**5, 2)
    L.append(ps_test1_2)
    ps_test1_3 = AKCN_ParameterSet(768, 2**11, 2**9, 2**4, 2)
    L.append(ps_test1_3)

    ps_test2_1 = AKCN_ParameterSet(1152, 2**12, 2**10, 2**3, 2)
    L.append(ps_test2_1)
    ps_test2_2 = AKCN_ParameterSet(1152, 2**13, 2**10, 2**3, 2)
    L.append(ps_test2_2)
 
    for ps_temp in L:
        print("**************************************")
        print("AKCN_RLWR: n = %d, q = 2^%d, p = 2^%d, g = 2^%d, eta = %d\n"%(ps_temp.n, log(ps_temp.q,2), log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta))  
        bandwidth(ps_temp)
        summarize(ps_temp)
        print()
    

    # ParameterSet(n, q, p, g, eta)
    # |K| = 40

    L_K320 = []

    #for k in range(,):
    #    for t in range(,):
    #        ps_test = AKCN_ParameterSet(, , , , , )
    #        L_K320.append(ps_test)

    ps_test1_K320 = AKCN_ParameterSet(768, 2**12, 2**9, 2**5, 2)
    L_K320.append(ps_test1_K320)
    ps_test2_K320 = AKCN_ParameterSet(768, 2**11, 2**9, 2**5, 2)
    L_K320.append(ps_test2_K320)
    ps_test3_K320 = AKCN_ParameterSet(768, 2**11, 2**9, 2**4, 2)
    L_K320.append(ps_test3_K320)

    for ps_temp in L_K320:
        print("**************************************")
        print("AKCN_RLWR: n = %d, q = 2^%d, p = 2^%d, g = 2^%d, eta = %d\n"%(ps_temp.n, log(ps_temp.q,2), log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta))
        bandwidth_K320(ps_temp)
        summarize(ps_temp)
        print()


    # ParameterSet(n, q, p, g, eta)
    # |K| = 64

    L_K512 = []

    #for k in range(,):
    #    for t in range(,):
    #        ps_test = AKCN_ParameterSet(, , , , , )
    #        L_K512.append(ps_test)
           
    ps_test1_K512 = AKCN_ParameterSet(1152, 2**12, 2**10, 2**3, 2)
    L_K512.append(ps_test1_K512)
    ps_test2_K512 = AKCN_ParameterSet(1152, 2**13, 2**10, 2**3, 2)
    L_K512.append(ps_test2_K512)

    for ps_temp in L_K512:
        print("**************************************")
        print("AKCN_RLWR: n = %d, q = 2^%d, p = 2^%d, g = 2^%d, eta = %d\n"%(ps_temp.n, log(ps_temp.q,2), log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta)) 
        bandwidth_K512(ps_temp)
        summarize(ps_temp)
        print()