
import math

from ld_estimator.utils import Haps

LN10 = math.log(10.0)

def likelihood_freqs(known, freqs, unknown):
    ''' get loglikelihood for model from estimated haplotype frequencies
    '''

    aa = freqs.aa
    ab = freqs.ab
    ba = freqs.ba
    bb = freqs.bb

    return (known.aa * math.log(aa) + \
        known.ab * math.log(ab) + \
        known.ba * math.log(ba) + \
        known.bb * math.log(bb) + \
        unknown * math.log(aa * bb + ab * ba)) / LN10

def likelihood_null(known, a1, a2, b1, b2, unknown):
    ''' get loglikelihood for null model from allele frequencies
    '''
    return (known.aa * math.log(a1 * a2) + \
        known.ab * math.log(a1 * b2) + \
        known.ba * math.log(b1 * a2) + \
        known.bb * math.log(b1 * b2) + \
        unknown * math.log(2 * a1 * a2 * b1 * b2)) / LN10

def get_lsurface(pA1, pA2, pB1, denom, known, unknown):
    ''' get log likelihood surface
    '''
    lsurface = []
    for i in range(100):
        dpr = i * 0.01
        AA = dpr * denom + pA1 * pA2
        AB = pA1 - AA
        BA = pA2 - AA
        BB = pB1 - BA
        freqs = Haps(AA, AB, BA, BB)

        if i == 100:
            # one value will be 0
            freqs.aa = max(1e-10, freqs.aa)
            freqs.ab = max(1e-10, freqs.ab)
            freqs.ba = max(1e-10, freqs.ba)
            freqs.bb = max(1e-10, freqs.bb)

        likelihood = likelihood_freqs(known, freqs, unknown)
        lsurface.append(likelihood)

    return lsurface
