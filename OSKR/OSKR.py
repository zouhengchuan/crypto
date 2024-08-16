from math import log, ceil
from OSKR_failure import p2_cyclotomic_error_probability, p2_cyclotomic_error_probability_AKCN, error_rate_test
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet
from proba_util import build_mod_switching_error_law

class OSKRParameterSet:
    def __init__(self, n, m, ks, ke,  q, rqk, rqc, rq2, ke_ct=None):
        if ke_ct is None:
            ke_ct = ke
        self.n = n
        self.m = m
        self.ks = ks     # binary distribution for the secret key
        self.ke = ke    # binary distribution for the ciphertext errors
        self.ke_ct = ke_ct    # binary distribution for the ciphertext errors
        self.q = q
        self.rqk = rqk  # 2^(bits in the public key)
        self.rqc = rqc  # 2^(bits in the first ciphertext)
        self.rq2 = rq2  # 2^(bits in the second ciphertext)


def OSKR_to_MLWE(kps):
    if kps.ks != kps.ke:
        raise "The security script does not handle different error parameter in secrets and errors (ks != ke) "

    # Check whether ciphertext error variance after rounding is larger than secret key error variance
    Rc = build_mod_switching_error_law(kps.q, kps.rqc)
    var_rounding = sum([i*i*Rc[i] for i in Rc.keys()])

    if kps.ke_ct/2. + var_rounding < kps.ke/2.:
        raise "The security of the ciphertext MLWE may not be stronger than the one of the public key MLWE"    

    return MLWEParameterSet(kps.n, kps.m, kps.m + 1, kps.ks, kps.q)


def communication_costs(ps):
    """ Compute the communication cost of a parameter set
    :param ps: Parameter set (ParameterSet)
    :returns: (cost_Alice, cost_Bob) (in Bytes)
    """
    A_space = ps.n/8. + ps.m * ceil(ps.n * log(ps.rqk)/log(2)/8.)   #pk
    B_space = ps.m * ceil(ps.n * log(ps.rqc)/log(2)/8.) + ceil(ps.n * log(ps.rq2)/log(2)/8.)  #ct
    return (int(A_space), int(B_space))


def summarize(ps):
    print ("params: ", ps.__dict__)
    print ("com costs: ", communication_costs(ps))

    h = error_rate_test(ps)
    F, f = p2_cyclotomic_error_probability(ps)
    print ("failure origin: %.1f = 2^%.5f"%(f, log(f + 2.**(-300))/log(2)))
    G, g = p2_cyclotomic_error_probability_AKCN(ps)
    print ("failure AKCN: %.1f = 2^%.5f"%(g, log(g + 2.**(-300))/log(2)))
    


if __name__ == "__main__":
    # Parameter sets
    ps_512 = OSKRParameterSet(256, 2, 3, 3, 3329, 2**12, 2**10, 2**4, ke_ct=2)
    ps_768 = OSKRParameterSet(256, 3, 2, 2, 3329, 2**12, 2**10, 2**4)
    ps_1024 = OSKRParameterSet(512, 2, 2, 2, 3329, 2**12, 2**11, 2**5)

    # Analyses
    print ("OSKR512 (light):")
    print ("--------------------")
    #print ("security:")
    #MLWE_summarize_attacks(Kyber_to_MLWE(ps_light))
    summarize(ps_512)
    print ()

    print ("OSKR768 (recommended):")
    print ("--------------------")
    #print ("security:")
    #MLWE_summarize_attacks(Kyber_to_MLWE(ps_recommended))
    summarize(ps_768)
    print ()

    print ("OSKR1024 (paranoid):")
    print ("--------------------")
    #print ("security:")
    #MLWE_summarize_attacks(Kyber_to_MLWE(ps_paranoid))
    summarize(ps_1024)
    print ()