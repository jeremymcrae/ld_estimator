
def get_allele_freqs(known, unknown):
    ''' calculate allele frequencies
    '''
    chroms = (known.aa + known.ab + known.ba + known.bb + (2 * unknown))
    pA1 = (known.aa + known.ab + unknown) / chroms
    pA2 = (known.aa + known.ba + unknown) / chroms

    return pA1, 1.0 - pA1, pA2, 1.0 - pA2
