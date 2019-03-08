
from ld_estimator.dprime import get_d

def fliphap(hap):
    ''' flip frequencies in a Haps object
    '''
    hap.aa, hap.ab = hap.ab, hap.aa
    hap.ba, hap.bb = hap.bb, hap.ba
    return hap

def flip_freqs(freqs, known, pA2, pB2):
    ''' flip frequencies to get positive D'
    '''
    # flip AA with AB and BA with BB
    freqs = fliphap(freqs)

    # flip frequency of second allele
    pA2, pB2 = pB2, pA2

    return freqs, fliphap(known), pA2, pB2, get_d(freqs)
