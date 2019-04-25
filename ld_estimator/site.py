
import pysam

from ld_estimator import pairwise_ld

def site_ld(vcf, chrom, pos, window=100000, subset=None):
    ''' calculate LD at one site to all surrounding variants

    Args:
        vcf: path to vcf, or pysam.VariantFile object
        chrom: chromosome.
        pos: start of region to calculate LD within
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

    variants = list(x for x in vcf.fetch(chrom, pos-1, pos) if x.pos == pos)
    if len(variants) == 0:
        raise ValueError(f'no variant found at: {chrom}:{pos}')

    v1 = variants[0]
    v1_geno = [v1.samples[x].alleles for x in v1.samples]

    lds = []
    for v2 in vcf.fetch(chrom, pos - window, pos + window):
        if v1 == v2:
            continue

        v2_geno = [v2.samples[x].alleles for x in v2.samples]
        ploidy = [False] * len(v2_geno)

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
