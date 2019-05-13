''' some common functions
'''

import pysam

from ld_estimator import pairwise_ld

PSEUDOAUTOSOMAL = {'grch37': {'X': [(60001, 2699520), (154931044, 155260560)],
                              'Y': [(10001, 2649520), (59034050, 59363566)]},
                   'grch38': {'X': [(10001, 2781479), (155701383, 156030895)],
                              'Y': [(10001, 2781479), (56887903, 57217415)]}}
MALE_CODES = set(['m', 'male', '1'])
FEMALE_CODES = set(['f', 'female', '2'])

def prepare_vcf(vcf, chrom, pos, subset, sexes, build='grch37'):
    ''' get a pysam.VariantFile ready for calculating LD and the ploidy state

    Opens the VCF file, and selects only the required samples (all samples by
    default).

    Args:
        vcf: path to VCF (or pysam.VariantFile handle)
        chrom: chromosome of interest
        pos: nucleotide position to check ploidy at
        subset: list of sample IDs to subset on (or None)
        sexes: list of sample sexes (or None). Has to match the subset order if
            that is present, otherwise the VCF samples order
        build: genome build to use (defaults to 'grch37')

    Returns:
        tuple of (pysam.VariantFile, ploidy list)
    '''
    try:
        vcf = pysam.VariantFile(vcf)
    except AttributeError:
        pass

    if not is_autosomal(chrom, pos) and sexes is None:
        raise ValueError('Sample sexes required for non-autosomal sites: {chrom}:{pos}')

    if sexes is None and subset is None:
        sexes = ['female'] * len(vcf.header.samples)

    if subset is not None:
        all_ids = set(vcf.header.samples)
        masked = [x in all_ids for x in subset]
        sample_ids = [sample for x, sample in zip(masked, subset) if x]
        sexes = [sample for x, sample in zip(masked, sexes) if x]
        try:
            vcf.subset_samples(sample_ids)
        except ValueError:
            # can only subset once, check the sets match if done previously
            assert all_ids == set(sample_ids)
    else:
        sample_ids = list(vcf.header.samples)

    # re-sort the sexes based on the VCF sample order, so that the ploidy list
    # later matches the genotypes order
    indexes = dict(zip(sample_ids, range(len(subset))))
    sample_ids = [sample_ids[indexes[x]] for x in vcf.header.samples]
    sexes = [sexes[indexes[x]] for x in vcf.header.samples]

    return vcf, get_ploidy(chrom, pos, sample_ids, sexes, build)

def check_ld(var1, var1_geno, var2, var2_geno, ploidy):
    ''' get a dictionary with variant coords and LD details
    '''
    data = {'chrom': var1.chrom,
            'var1': {'pos': var1.pos, 'id': var1.id, 'ref': var1.ref, 'alt': var1.alts[0]},
            'var2': {'pos': var2.pos, 'id': var2.id, 'ref': var2.ref, 'alt': var2.alts[0]},
            'r_squared': None, 'dprime': None, 'phase': None}

    linkage = pairwise_ld(var1_geno, var2_geno, ploidy)
    if linkage is not None:
        data['r_squared'] = linkage.r_squared
        data['dprime'] = linkage.dprime
        data['phase'] = linkage.phase

    return data

def is_autosomal(chrom, pos, build='grch37', par_regions=None):
    ''' identify whether a genomic site is autosomal

    This includes pseudoautosomal regions as autosomal.

    Args:
        chrom: chromosome of variant
        pos: nucleotide position of variant
        build: genome build (defaults to 'grch37', but can also use 'grch38',
            'hg19', 'hg38' and use the default pseudoautosomal regions. If a
            a different genome build is needed, supply a own par_regions
            dictionary.
        par_regions: dict of pseudoautosomal regions, indexed by genome build,
            then chromosome. Regions are given as list of (start, end) tuples.
            For example {'grch37': {'X': [(1, 10), ...], 'Y': [(30, 40), ...]}}
    '''
    # standardise the genome build name
    build = build.lower()
    names = {'hg19': 'grch37', 'hg38': 'grch38'}
    build = names[build] if build in names else build

    # load the default PAR dict if none provided
    if not par_regions:
        par_regions = PSEUDOAUTOSOMAL

    # select the pseudoautosomal regions on the given chromosome
    chrom = chrom.strip('chr')
    par = par_regions[build]
    autosomal = chrom not in par or any(x[0] <= int(pos) <= x[1] for x in par[chrom])
    return autosomal

def get_ploidy(chrom, pos, samples, sexes, build='grch37'):
    ''' figure out the expected ploidy for each individual

    Given a list of sexes, and the chromosome and position, figure out whether
    the individuals are haploid (males on chrX). This ignores chrY, but that is
    ok, since chrY is haploid, and therefore LD is simple to assess.

    Currently this only checks ploidy at one site, but technically we should
    figure out ploidy at both sites being assessed for LD, but that gets tricky.
    '''
    if is_autosomal(chrom, pos, build):
        return [False] * len(samples)

    for sample, sex in zip(samples, sexes):
        sex = sex.lower()
        assert sex in MALE_CODES or sex in FEMALE_CODES, f'unknown sex ({sex}) in {sample}'

    return [sex.lower() in MALE_CODES for _, sex in sexes]
