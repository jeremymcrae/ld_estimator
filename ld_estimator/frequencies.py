
from ld_estimator.utils import Haps

def count_haplotypes(iteration: int, known: Haps, freqs: Haps, unknown: int) -> Haps:
    ''' count haplotypes
    '''
    counts = Haps(known.aa, known.ab, known.ba, known.bb)
    if iteration > 0:
        obligates = freqs.aa * freqs.bb
        hets = freqs.ab * freqs.ba
        counts.aa += unknown * obligates / (obligates + hets)
        counts.bb += unknown * obligates / (obligates + hets)
        counts.ab += unknown * hets / (obligates + hets)
        counts.ba += unknown * hets / (obligates + hets)

    return counts

def estimate_frequencies(counts: Haps, prob: float) -> Haps:
    ''' estimate haplotype frequencies from haplotype counts
    '''
    total = sum(counts) + (4.0 * prob)

    # estimate frequency of each haplotype, but keep them above 1e-10
    aa = max(1e-10, ((counts.aa + prob) / total))
    ab = max(1e-10, ((counts.ab + prob) / total))
    ba = max(1e-10, ((counts.ba + prob) / total))
    bb = max(1e-10, ((counts.bb + prob) / total))

    return Haps(aa, ab, ba, bb)
