''' provides function to calculate LD within a genome region, given a VCF
'''

from itertools import combinations

from ld_estimator.utils import prepare_vcf, check_ld
from ld_estimator.pairwise import _is_monomorphic

def region_ld(vcf, chrom, start, end, subset=None, sexes=None, build='grch37'):
    ''' calculate LD within a genome region

    Args:
        vcf: path to vcf, or pysam.VariantFile object
        chrom: chromosome.
        start: start of region to calculate LD within
        end: end of region to calculate LD within
        subset: list of sample IDs to include in LD calculation. Default is to
            include every sample.
        sexes: list of sexes of samples (matching subset order). Only needed
            if variant on a sex chromosome and not in a pseudoautosomal region.
        build: genome build to use (defaults to 'grch37')

    Returns:
        list of LD estimates, each as a dict with [chrom, var1_pos, var2_pos,
            r_squared and dprime] keys.
    '''

    vcf, ploidy = prepare_vcf(vcf, chrom, start, subset, sexes, build)

    lds = []

    variants = []
    for var in vcf.fetch(chrom, start, end):
        geno = [var.samples[x].alleles for x in var.samples]
        if _is_monomorphic(geno):
            continue
        geno = [tuple(map(str.encode, x)) for x in geno]
        variants.append((var, geno))

    for (var1, var1_geno), (var2, var2_geno) in combinations(variants, 2):
        data = check_ld(var1, var1_geno, var2, var2_geno, ploidy)
        lds.append(data)

    return lds
