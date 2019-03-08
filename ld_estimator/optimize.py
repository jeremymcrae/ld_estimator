
from ld_estimator.utils import Haps
from ld_estimator.frequencies import count_haplotypes, estimate_frequencies
from ld_estimator.likelihoods import likelihood_freqs

def get_frequencies(known, unknown, epsilon=0.00000001):
    ''' run expectation maximisation to get optimal haplotype frequencies
    '''
    initial = Haps(0.01, 0.01, 0.01, 0.01)
    counts = count_haplotypes(0, known, initial, unknown)
    frequencies = estimate_frequencies(counts, 0.1)

    iteration = 1
    current = -999999999.0
    previous = 5
    while abs(current - previous) > epsilon and iteration < 1000:
        previous = current
        counts = count_haplotypes(iteration, known, frequencies, unknown)
        current = likelihood_freqs(known, frequencies, unknown)
        frequencies = estimate_frequencies(counts, 0.0)
        iteration += 1

    return frequencies
