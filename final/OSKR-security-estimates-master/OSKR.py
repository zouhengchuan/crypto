from math import log, ceil
from OSKR_failure import*
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from proba_util import build_mod_switching_error_law

class OSKRParameterSet:
    def __init__(self, n, m, ks, ke,  q, rqc, rq2, ke_ct=None):
        if ke_ct is None:
            ke_ct = ke
        self.n = n
        self.m = m
        self.ks = ks        # binary distribution for the secret key
        self.ke = ke        # binary distribution for the errors
        self.ke_ct = ke_ct  # binary distribution for the ciphertext errors
        self.q = q
        self.rqc = rqc      # 2^(bits in the first ciphertext)
        self.rq2 = rq2      # 2^(bits in the second ciphertext)
        self.sigma = sqrt(sqrt(ks * ke )/2)


def OSKR_to_MLWE(kps):
    Rc = build_mod_switching_error_law(kps.q, kps.rqc)
    var_rounding = sum([i*i*Rc[i] for i in Rc.keys()])

    if kps.ke_ct/2. + var_rounding < kps.ke/2.:
        raise "The security of the ciphertext MLWE may not be stronger than the one of the public key MLWE"    
    return MLWEParameterSet(kps.n * kps.m, kps.q, kps.sigma)


def summarize(ps):
    print ("params: ", ps.__dict__)
    if ps.n == 384:
        bandwidth_K320(ps)
    else:
        bandwidth_K512(ps)
    error_probability(ps)
    error_probability_1(ps)
    


if __name__ == "__main__":

    L = []
    eta_k = 3
    eta_e = 3
    eta_e_ct = 2
    for u in range(11,10,-1):
        for v in range(7,3,-1):
            ps_test = OSKRParameterSet(576, 2, eta_k, eta_e, 6337, 2**u, 2**v, ke_ct = eta_e_ct)
            L.append(ps_test)

    # ps = OSKRParameterSet(576, 2, 3, 3, 6337, 2**11, 2**7, ke_ct=2)
    # ps2 = OSKRParameterSet(384, 2, 3, 3, 6337, 2**11, 2**3, ke_ct=2)
    # L = [ps,ps2]
    for ps_test in L:
        print ("--------------------")
        print ("security:")
        MLWE_summarize_attacks(OSKR_to_MLWE(ps_test))
        summarize(ps_test)
        print ()