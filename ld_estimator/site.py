
import pysam

from ld_estimator import pairwise_ld

def site_ld(vcf, chrom, pos, window=100000, subset=None):
    ''' calculate LD at one site to all surrounding variants

    Args:
        vcf: path to vcf, or pysam.VariantFile object
        chrom: chromosome.
        start: start of region to calculate LD within
        end: end of region to calculate LD within
        window: number of base pairs to check LD within
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

    v1 = list(vcf.fetch(chrom, pos-1, pos))
    v1_geno = [v1.samples[x].alleles for x in v1.samples]

    lds = []
    for v2 in vcf.fetch(chrom, pos - window, pos + window):
        if v1 == v2:
            continue

        v2_geno = [v2.samples[x].alleles for x in v2.samples]
        ploidy = [False] * len(v2_geno)

        data = {'chrom': chrom, 'var1_pos': v1.pos, 'var1_id': v1.id,
            'var2_pos': v2.pos, 'var2_id': v2.id}

        ld = pairwise_ld(v1_geno, v2_geno, ploidy)
        if ld is None:
            data['r_squared'] = None
            data['dprime'] = None
        else:
            data['r_squared'] = ld.r_squared
            data['dprime'] = ld.dprime

        lds.append(data)

    return lds
