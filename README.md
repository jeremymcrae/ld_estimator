## ld_estimator: A package to estimate linkage disequilibrium
Computes linkage disequilibrium in python. This uses the maximum-likelihood of
[Excoffier & Slatkin](https://doi.org/10.1093/oxfordjournals.molbev.a040269),
just converted from [HaploView](https://github.com/jazzywhit/Haploview) to c++
, with python bindings. It's not too slow, can calculate ~1000 pairs per
second.

### Installation
The simplest way to install ld_estimator is through pip:
```sh
pip install ld_estimator
```

### Usage
Use ld_estimator within a python environment
```python
from ld_estimator import pairwise_ld

var1 = [(0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (0, 1), (0, 0), (0, 0), (1, 1)]
var2 = [(0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (0, 1), (0, 0), (1, 1), (1, 1)]
is_haploid = [False, False, False, False, False, False, False, True, True, True]
ld = pairwise_ld(var1, var2, is_haploid)
print(ld.dprime)
print(ld.r_squared)

# or calculate LD for all pairs of variants in a region in a VCF:
from ld_estimator import region_ld
vcf_path = 'PATH_TO_VCF'
ld = region_ld(vcf_path, chrom, start, end)

# or calculate LD to a site within a region in a VCF. This defaults to checking
# variants within a 100 kb window of the specified site.
from ld_estimator import site_ld
vcf_path = 'PATH_TO_VCF'
ld = site_ld(vcf_path, chrom, pos, window=200000)

# can pass in multiple positions in the same region at once
ld = site_ld(vcf_path, chrom, [pos2, pos2, pos3], window=200000)

# both region_ld() and site_ld() can take a list of sample IDs to subset the
# samples used for calculating LD. For example:
ld = site_ld(vcf_path, chrom, pos, subset=['sample1', 'sample2'])

# if the variant is on a sex chromosome, you'll have to pass in a list of sample
# sexes (matching order of the subset IDs if present, otherwise the VCF samples)
ld = site_ld(vcf_path, 'X', 20000000, sexes=['male', 'female'])
```
