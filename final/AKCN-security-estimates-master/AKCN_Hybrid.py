from math import sqrt, log
from proba_util import *
from AKCN_Hybrid_failure import *
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet


class AKCN_ParameterSet:
    def __init__(self, n, q, p, g, eta1, eta2, eta):
        self.n = n
        self.q = q
        self.p = p
        self.g = g
        self.eta1 = eta1
        self.eta2 = eta2
        self.eta = eta
        self.sigmaE1 = sqrt(var_of_law(build_centered_binomial_law(eta1)))
        self.sigmaE2 = sqrt(var_of_law(build_centered_binomial_law(eta2)))
        self.sigmaR1 = sqrt(var_of_law(build_centered_binomial_law(eta)))
        self.sigmaR2 = sqrt(var_of_law(build_rounding_law_rlwr(q,p)))

        # standard deviation used for RLWE
        self.sigmaE = sqrt(self.sigmaE1 * self.sigmaE2)

        # standard deviation used for RLWR
        self.sigmaR = sqrt(self.sigmaR1 * self.sigmaR2)
        
        

def AKCN_to_MLWE(ps):
    return MLWEParameterSet(ps.n, ps.q, ps.sigmaE)

def AKCN_to_MLWE_RLWR(ps):
    return MLWEParameterSet(ps.n, ps.q, ps.sigmaR)

def summarize(ps):
    print("Security of RLWE")
    MLWE_summarize_attacks(AKCN_to_MLWE(ps))
    print("\nSecurity of RLWR")
    MLWE_summarize_attacks(AKCN_to_MLWE_RLWR(ps))
    ErrorRate_Hybrid(ps)
    

if __name__ == "__main__":

    # ParameterSet(n, q, p, g, eta1, eta2, eta)
    # |K| = n/16

    L = []

    ps_test1_1 = AKCN_ParameterSet(768, 2**11, 2**9, 2**4, 2, 2, 2)
    L.append(ps_test1_1)
    ps_test2_1 = AKCN_ParameterSet(1152, 2**12, 2**10, 2**3, 2, 2, 2)
    L.append(ps_test2_1)
 
    for ps_temp in L:
        print("**************************************")
        print("AKCN_Hybrid: n = %d, q = %d, p = 2^%d, g = 2^%d, eta1 = %d, eta2 = %d, eta = %d\n"%(ps_temp.n, ps_temp.q, log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta1, ps_temp.eta2, ps_temp.eta))  
        bandwidth(ps_temp)
        summarize(ps_temp)
        print()
    


    # ParameterSet(n, q, p, g, eta1, eta2, eta)
    # |K| = 40

    L_K320 = []

    #for k in range(,):
    #    for t in range(,):
    #        ps_test = AKCN_ParameterSet(, , , , , )
    #        L_K320.append(ps_test)

    ps_test1_K320 = AKCN_ParameterSet(768, 2**11, 2**9, 2**4, 2, 2, 2)
    L_K320.append(ps_test1_K320)

    for ps_temp in L_K320:
        print("**************************************")
        print("AKCN_Hybrid: n = %d, q = %d, p = 2^%d, g = 2^%d, eta1 = %d, eta2 = %d, eta = %d\n"%(ps_temp.n, ps_temp.q, log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta1, ps_temp.eta2, ps_temp.eta))  
        bandwidth_K320(ps_temp)
        summarize(ps_temp)
        print()


    # ParameterSet(n, q, p, g, eta1, eta2, eta)
    # |K| = 64

    L_K512 = []

    #for k in range(,):
    #    for t in range(,):
    #        ps_test = AKCN_ParameterSet(, , , , , )
    #        L_K512.append(ps_test)

    ps_test1_K512 = AKCN_ParameterSet(1152, 2**12, 2**10, 2**3, 2, 2, 2)
    L_K512.append(ps_test1_K512)

    for ps_temp in L_K512:
        print("**************************************")
        print("AKCN_Hybrid: n = %d, q = %d, p = 2^%d, g = 2^%d, eta1 = %d, eta2 = %d, eta = %d\n"%(ps_temp.n, ps_temp.q, log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta1, ps_temp.eta2, ps_temp.eta))  
        bandwidth_K512(ps_temp)
        summarize(ps_temp)
        print()
    
