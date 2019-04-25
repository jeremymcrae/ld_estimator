
from itertools import combinations

import pysam

from ld_estimator import pairwise_ld
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

    try:
        vcf = pysam.VariantFile(vcf)
    except AttributeError:
        pass

    if subset is not None:
        subset = [x for x in subset if x in vcf.header.samples]
        vcf.subset_samples(subset)

    lds = []

    vars = []
    for var in vcf.fetch(chrom, start, end):
        geno = [var.samples[x].alleles for x in var.samples]
        if _is_monomorphic(geno):
            continue
        vars.append((var, geno))

    ploidy = [False] * len(vars[0][1])
    for (v1, v1_geno), (v2, v2_geno) in combinations(vars, 2):
        data = {'chrom': chrom, 'var1_pos': v1.pos, 'var1_id': v1.id,
            'var2_pos': v2.pos, 'var2_id': v2.id, 'r_squared': None,
            'dprime': None, 'phase': None}

        ld = pairwise_ld(v1_geno, v2_geno, ploidy)
        if ld is not None:
            data['r_squared'] = ld.r_squared
            data['dprime'] = ld.dprime
            data['phase'] = ld.phase

        lds.append(data)

    return lds
