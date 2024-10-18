from math import sqrt, log
from proba_util import *
from AKCN_RLWE_failure import *
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
        self.sigma1 = sqrt(eta1/2.)
        self.sigma2 = sqrt(eta2/2.)
        self.sigma = sqrt(eta/2.)

        # standard deviation used in LWE attack
        self.sigmaE = sqrt(sqrt(eta1*eta2)/2)

        

def AKCN_to_MLWE(ps):

    # Check whether ciphertext error variance after rounding is larger than secret key error variance
    Rc = build_mod_switching_error_law(ps.q, ps.p)
    var_rounding = sum([i*i*Rc[i] for i in Rc.keys()])

    if ps.eta/2. + var_rounding < ps.eta2/2.:
        raise "The security of the ciphertext MLWE may not be stronger than the one of the public key MLWE"    

    return MLWEParameterSet(ps.n, ps.q, ps.sigmaE) 


def summarize(ps):
    MLWE_summarize_attacks(AKCN_to_MLWE(ps))
    ErrorRate_RLWE(ps)


if __name__ == "__main__":


    # ParameterSet(n, q, p, g, eta_1, eta_2, eta)
    # |K| = n/16

    L = []

    #for k in range(,):
    #    for t in range(,):
    #       ps_test = AKCN_ParameterSet(, , , , , , )
    #       L.append(ps_test)

    # ps_test1_1 = AKCN_ParameterSet(768, 3457, 2**10, 2**3, 2, 2, 2)
    # L.append(ps_test1_1)
    # ps_test1_2 = AKCN_ParameterSet(768, 3457, 2**10, 2**3, 2, 3, 3)
    # L.append(ps_test1_2)
    # ps_test2_1 = AKCN_ParameterSet(1152, 3457, 2**10, 2**4, 2, 2, 2)
    # L.append(ps_test2_1)
    # ps_test2_2 = AKCN_ParameterSet(1152, 3457, 2**10, 2**5, 2, 3, 3)
    # L.append(ps_test2_2)

    # for ps_temp in L:
    #     print("**************************************")
    #     print("AKCN_RLWE: n = %d, q = %d, p = 2^%d, g = 2^%d, eta1 = %d, eta2 = %d, eta = %d\n"%(ps_temp.n, ps_temp.q, log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta1, ps_temp.eta2, ps_temp.eta))  
    #     bandwidth(ps_temp)
    #     summarize(ps_temp)
    #     print()
    

    # ParameterSet(n, q, p, g, eta_1, eta_2, eta)
    # |K| = 40

    L_K320 = []

    #for k in range(,):
    #    for t in range(,):
    #       ps_test = AKCN_ParameterSet(, , , , , , )
    #       L_K320.append(ps_test)

    # ps_test1_K320 = AKCN_ParameterSet(768, 3457, 2**10, 2**3, 2, 2, 2)
    # L_K320.append(ps_test1_K320)
    ps_test2_K320 = AKCN_ParameterSet(768, 3457, 2**10, 2**3, 2, 3, 2)
    L_K320.append(ps_test2_K320)

    for ps_temp in L_K320:
        print("**************************************")
        print("AKCN_RLWE: n = %d, q = %d, p = 2^%d, g = 2^%d, eta1 = %d, eta2 = %d, eta = %d\n"%(ps_temp.n, ps_temp.q, log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta1, ps_temp.eta2, ps_temp.eta))  
        bandwidth_K320(ps_temp)
        summarize(ps_temp)
        print()


    # ParameterSet(n, q, p, g, eta_1, eta_2, eta)
    # |K| = 64

    L_K512 = []

    #for k in range(,):
    #    for t in range(,):
    #       ps_test = AKCN_ParameterSet(, , , , , , )
    #       L_K512.append(ps_test)

    # ps_test1_K512 = AKCN_ParameterSet(1152, 3457, 2**10, 2**4, 2, 2, 2)
    # L_K512.append(ps_test1_K512)
    ps_test2_K512 = AKCN_ParameterSet(1152, 3457, 2**10, 2**5, 2, 3, 2)
    L_K512.append(ps_test2_K512)

    for ps_temp in L_K512:
        print("**************************************")
        print("AKCN_RLWE: n = %d, q = %d, p = 2^%d, g = 2^%d, eta1 = %d, eta2 = %d, eta = %d\n"%(ps_temp.n, ps_temp.q, log(ps_temp.p,2), log(ps_temp.g,2), ps_temp.eta1, ps_temp.eta2, ps_temp.eta))  
        bandwidth_K512(ps_temp)
        summarize(ps_temp)
        print()