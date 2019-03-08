
from collections import Counter

class Haps(object):
    __slots__ = ['aa', 'ab', 'ba', 'bb']
    def __init__(self, aa, ab, ba, bb):
        self.aa = aa
        self.ab = ab
        self.ba = ba
        self.bb = bb
    def __repr__(self):
        return f'Haps({self.aa}, {self.ab}, {self.ba}, {self.bb})'
    def __iter__(self):
        for x in [self.aa, self.ab, self.ba, self.bb]:
            yield x

def is_monomorphic(var):
    ''' check that a variant has two or more alleles
    '''
    nones = set([None])
    alleles = set()
    for geno in var:
        alleles |= (set(geno) - nones)
        if len(alleles) > 1:
            return False

    return True

def lacks_haplotypes(counts, unknown):
    ''' check if the two variants lack some haplotypes
    '''
    r1 = counts.aa + counts.ab
    r2 = counts.ba + counts.bb
    c1 = counts.aa + counts.ba
    c2 = counts.ab + counts.bb
    return (r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0) and unknown == 0

def get_alleles(var):
    ''' get the major and minor alleles for a variant
    
    The minor allele is the second most common allele in the population.
    We've previously excluded non-polymorphic sites, so the variant has at
    least two alleles.
    '''
    counts = Counter(allele for geno in var for allele in geno)

    # remove missing alleles
    if None in counts:
        del counts[None]

    alleles = sorted(counts, key=lambda x: counts[x], reverse=True)

    return alleles[0], alleles[1]
