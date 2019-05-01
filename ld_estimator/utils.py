''' some common functions
'''

import pysam

from ld_estimator import pairwise_ld

def prepare_vcf(vcf, subset):
    ''' get a pysam.VariantFile ready for calculating LD
    '''
    try:
        vcf = pysam.VariantFile(vcf)
    except AttributeError:
        pass

    if subset is not None:
        subset = [x for x in subset if x in vcf.header.samples]
        vcf.subset_samples(subset)

    return vcf


def check_ld(var1, var1_geno, var2, var2_geno, ploidy):
    ''' get a dictionary with variant coords and LD details
    '''
    data = {'chrom': var1.chrom, 'var1_pos': var1.pos, 'var1_id': var1.id,
            'var2_pos': var2.pos, 'var2_id': var2.id, 'r_squared': None,
            'dprime': None, 'phase': None}

    linkage = pairwise_ld(var1_geno, var2_geno, ploidy)
    if linkage is not None:
        data['r_squared'] = linkage.r_squared
        data['dprime'] = linkage.dprime
        data['phase'] = linkage.phase

    return data
