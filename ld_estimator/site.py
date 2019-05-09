''' provides function to calculate LD against a single site in a VCF
'''

from ld_estimator.utils import prepare_vcf, check_ld

def site_ld(vcf, chrom, positions, window=100000, subset=None, sexes=None, build='grch37'):
    ''' calculate LD at specified site/s to all surrounding variants

    Args:
        vcf: path to vcf, or pysam.VariantFile object
        chrom: chromosome.
        positions: position/s of variant/s to calculate LD against.
        window: number of base pairs to check LD within
        subset: list of sample IDs to include in LD calculation.
            Default is to include every sample.
        sexes: list of sexes of samples (matching subset order). Only needed
            if variant on a sex chromosome and not in a pseudoautosomal region.
        build: genome build to use (defaults to 'grch37')

    Returns:
        list of LD estimates, each as a dict with [chrom, var1_pos, var2_pos,
            var1_id, var2_id, r_squared, dprime and phase] keys.
    '''
    if not isinstance(positions, list):
        positions = [positions]

    vcf, ploidy = prepare_vcf(vcf, chrom, positions[0], subset, sexes, build)

    variants = []
    for pos in positions:
        entries = list(x for x in vcf.fetch(chrom, pos-1, pos) if x.pos == pos)
        if entries == []:
            raise ValueError(f'no variant found at: {chrom}:{pos}')
        var = entries[0]
        geno = [var.samples[x].alleles for x in var.samples]
        geno = [tuple(map(str.encode, x)) for x in geno]
        variants.append([var, geno])

    lds = []
    for var2 in vcf.fetch(chrom, min(positions) - window, max(positions) + window):
        var2_geno = [var2.samples[x].alleles for x in var2.samples]
        for var1, var1_geno in variants:
            if var1 == var2:
                continue
            if abs(var1.pos - var2.pos) > window:
                continue

            data = check_ld(var1, var1_geno, var2, var2_geno, ploidy)
            lds.append(data)

    return lds
