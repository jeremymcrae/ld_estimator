''' provides function to calculate LD against a single site in a VCF
'''

from ld_estimator.utils import prepare_vcf, check_ld

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

    vcf = prepare_vcf(vcf, subset)

    variants = list(x for x in vcf.fetch(chrom, pos-1, pos) if x.pos == pos)
    if variants == []:
        raise ValueError(f'no variant found at: {chrom}:{pos}')

    var1 = variants[0]
    var1_geno = [var1.samples[x].alleles for x in var1.samples]

    lds = []
    for var2 in vcf.fetch(chrom, pos - window, pos + window):
        if var1 == var2:
            continue

        var2_geno = [var2.samples[x].alleles for x in var2.samples]
        ploidy = [False] * len(var2_geno)
        data = check_ld(var1, var1_geno, var2, var2_geno, ploidy)
        lds.append(data)

    return lds
