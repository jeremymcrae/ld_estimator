
from collections import defaultdict

from ld_estimator.utils import get_alleles, Haps

def count_haplotypes(var1, var2, ploidy):
    ''' tally known haplotypes of known and unknown phase

    Args:
        var1: list of genotypes for first variant
        var2: list of genotypes for second variant, ordered as per var1
        ploidy: list of ploidy states, ordered as per var1
    '''
    unknown = 0
    major_1, minor_1 = get_alleles(var1)
    major_2, minor_2 = get_alleles(var2)

    slots = defaultdict(lambda: defaultdict(int))
    slots[major_1][major_2] = 0
    slots[major_1][minor_2] = 0
    slots[minor_1][major_2] = 0
    slots[minor_1][minor_2] = 0

    # iterate through all chromosomes in dataset
    for (a1, a2), (b1, b2), is_haploid in zip(var1, var2, ploidy):
        if not is_haploid:
            if None in [a1, a2, b1, b2]:
                # skip missing data
                continue
            elif a1 != a2 and b1 != b2:
                # both genotypes are heterozygous
                unknown += 1
            elif a1 != a2:
                # b alleles are homozygous
                slots[major_1][b1] += 1
                slots[minor_1][b1] += 1
            elif b1 != b2:
                # a alleles are homozygous
                slots[a1][major_2] += 1
                slots[a1][minor_2] += 1
            else:
                # both are homozygous
                slots[a1][b1] += 2
        else:
            # include haploid chromosomes
            if a1 is not None and a2 is not None:
                slots[a1][b1] += 1

    aa = slots[major_1][major_2]
    ab = slots[major_1][minor_2]
    ba = slots[minor_1][major_2]
    bb = slots[minor_1][minor_2]
    counts = Haps(aa, ab, ba, bb)

    return counts, unknown
