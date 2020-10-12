# cython: language_level=3, boundscheck=False

from libc.stdlib cimport malloc, free
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

from cython.operator cimport dereference as deref

from ld_estimator.linkage import LD

cdef extern from 'haps.h' namespace 'ld_estimator':
  cdef cppclass Haps[T]:
      T aa
      T ab
      T ba
      T bb

cdef extern from 'linkage.h' namespace 'ld_estimator':
  cdef struct Phase:
      string var1_allele
      string var2_allele
  cdef struct Linkage:
      double dprime
      double loglikelihood
      double r_squared
      Haps[double] freqs
      Phase phase

cdef extern from 'ld.h' namespace 'ld_estimator':
    Linkage pairwise(char **a1, char **a2, char **b1, char **b2, int size, vector[bool] &) except +

cdef to_bytes(var):
    return [[str(a).encode('utf8') if a is not None else b'' for a in g] for g in var]

cdef to_cstring_array(list_str, char ** a1, char ** a2, idx1, idx2):
    for i, x in enumerate(list_str):
        a1[i] = x[idx1]
        a2[i] = x[idx2]

def pairwise_ld(var1, var2, ploidy):
    ''' check LD for pair of variants

    Args:
        var1: alleles for first variant, as list of (allele1, allele2) lists per sample
        var2: alleles for second variant, as list of (allele1, allele2) lists
            per sample. Sample order must be same as for var1
        ploidy: list of haploid sample states, same length as var1 and var2 and
            same sample order.

    Returns:
        LD obect, with dprime, loglikelihood, r-squared, [freqs] and [phases] attributes
    '''
    # only convert to allele vectors to bytes if not done already
    if not isinstance(var1[0][0], bytes):
        var1 = to_bytes(var1)
    if not isinstance(var2[0][0], bytes):
        var2 = to_bytes(var2)

    if len(var1) != len(var2):
        raise ValueError('genotypes lists need to be the same length')
    if len(var1) == 0:
        raise ValueError('zero genotypes in lists supplied for LD')

    cdef char **a1 = <char **>malloc(len(var1) * sizeof(char *))
    cdef char **a2 = <char **>malloc(len(var1) * sizeof(char *))
    cdef char **b1 = <char **>malloc(len(var1) * sizeof(char *))
    cdef char **b2 = <char **>malloc(len(var1) * sizeof(char *))

    to_cstring_array(var1, a1, a2, 0, -1)
    to_cstring_array(var2, b1, b2, 0, -1)
    try:
        ld = pairwise(a1, a2, b1, b2, len(var1), ploidy)
    except ValueError:
        return None

    free(a1)
    free(a2)
    free(b1)
    free(b2)

    return LD(ld.dprime, ld.loglikelihood, ld.r_squared,
      [ld.freqs.aa, ld.freqs.ab, ld.freqs.ba, ld.freqs.bb], [ld.phase.var1_allele.decode('utf8'), ld.phase.var2_allele.decode('utf8')])

cdef extern from 'utils.h' namespace 'ld_estimator':
    bool is_monomorphic(char **a1, char **a2, int size) except +

def _is_monomorphic(var):
    cdef char **a1 = <char **>malloc(len(var) * sizeof(char *))
    cdef char **a2 = <char **>malloc(len(var) * sizeof(char *))
    to_cstring_array(var, a1, a2, 0, -1)

    mono = is_monomorphic(a1, a2, len(var))

    free(a1)
    free(a2)

    return mono
