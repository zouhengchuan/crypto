from math import log
import operator as op
from math import factorial as fac
from math import sqrt, log, ceil
from proba_util import *
from AKCN_RLWR_failure import *
from MLWE_security import MLWE_summarize_attacks, MLWEParameterSet



sigma1 = sqrt(var_of_law(build_centered_binomial_law(3)))

sigma2 = sqrt(var_of_law(build_rounding_law_rlwr(12,10)))

sigma = sqrt(sigma1 * sigma2)



MLWE_summarize_attacks(MLWEParameterSet(1024, 2**12, sigma))