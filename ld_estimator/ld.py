
from ld_estimator.tallies import count_haplotypes
from ld_estimator.confidence_interval import ld_confidence_interval
from ld_estimator.likelihoods import likelihood_freqs, likelihood_null, get_lsurface
from ld_estimator.optimize import get_frequencies
from ld_estimator.dprime import get_d, get_denominator
from ld_estimator.flip import flip_freqs
from ld_estimator.linkage import LD
from ld_estimator.af import get_allele_freqs
from ld_estimator.utils import is_monomorphic, lacks_haplotypes

def estimate_ld(var1, var2, ploidy):
    ''' compute linkage disequilibrium between two variants

    Args:
        var1: list of genotypes for first variant
        var2: list of genotypes for second variant, ordered as per var1
        ploidy: list of ploidy states, ordered as per var1
    '''
    if is_monomorphic(var1) or is_monomorphic(var2):
        return None

    known, unknown = count_haplotypes(var1, var2, ploidy)
    if lacks_haplotypes(known, unknown):
        return LD(1, 0, 0, 0, 0, [0])

    pA1, pB1, pA2, pB2 = get_allele_freqs(known, unknown)

    freqs = get_frequencies(known, unknown)
    loglike1 = likelihood_freqs(known, freqs, unknown)
    loglike0 = likelihood_null(known, pA1, pA2, pB1, pB2, unknown)

    d = get_d(freqs)
    if d < 0:
        freqs, known, pA2, pB2, d = flip_freqs(freqs, known, pA2, pB2)

    denom = get_denominator(freqs)

    rsq = (d * d) / (pA1 * (1.0 - pA1) * pA2 * (1.0 - pA2))
    surface = get_lsurface(pA1, pA2, pB1, denom, known, unknown)
    low_ci, high_ci = ld_confidence_interval(surface, loglike1, alpha=0.05)

    return LD(d / denom, loglike1 - loglike0, rsq, low_ci, high_ci, list(freqs))
