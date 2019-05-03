''' provides function to calculate LD within a genome region, given a VCF
'''

from itertools import combinations

from ld_estimator.utils import prepare_vcf, check_ld
from ld_estimator.pairwise import _is_monomorphic

def region_ld(vcf, chrom, start, end, subset=None):
    ''' calculate LD within a genome region

    Args:
        vcf: path to vcf, or pysam.VariantFile object
        chrom: chromosome.
        start: start of region to calculate LD within
        end: end of region to calculate LD within
        subset: list of sample IDs in VCF to include in LD calculation. Default
            is to include every sample.

    Returns:
        list of LD estimates, each as a dict with [chrom, var1_pos, var2_pos,
            r_squared and dprime] keys.
    '''

    vcf = prepare_vcf(vcf, subset)

    lds = []

    variants = []
    for var in vcf.fetch(chrom, start, end):
        geno = [var.samples[x].alleles for x in var.samples]
        if _is_monomorphic(geno):
            continue
        geno = [tuple(map(str.encode, x)) for x in geno]
        variants.append((var, geno))

    ploidy = [False] * len(variants[0][1])
    for (var1, var1_geno), (var2, var2_geno) in combinations(variants, 2):
        data = check_ld(var1, var1_geno, var2, var2_geno, ploidy)
        lds.append(data)

    return lds
