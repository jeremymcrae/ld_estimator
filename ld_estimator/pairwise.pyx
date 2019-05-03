# cython: language_level=3, boundscheck=False

from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

from ld_estimator.linkage import LD

cdef extern from 'haps.h' namespace 'ld_estimator':
  cdef cppclass Haps[T]:
      T aa
      T ab
      T ba
      T bb

cdef extern from 'linkage.h' namespace 'ld_estimator':
  cdef struct Linkage:
      double dprime
      double loglikelihood
      double r_squared
      double ci_low
      double ci_high
      Haps[double] freqs
      vector[string] phase

cdef extern from 'ld.h' namespace 'ld_estimator':
    Linkage pairwise(vector[vector[string]], vector[vector[string]], vector[bool]) except +

cdef to_bytes(var):
    return [[str(a).encode('utf8') for a in g] for g in var]

def pairwise_ld(var1, var2, vector[bool] ploidy):
    # only convert to allele vectors to bytes if not done already
    if not isinstance(var1[0][0], bytes):
        var1 = to_bytes(var1)
    if not isinstance(var2[0][0], bytes):
        var2 = to_bytes(var2)
    
    try:
        ld = pairwise(var1, var2, ploidy)
    except ValueError:
        return None
    
    return LD(ld.dprime, ld.loglikelihood, ld.r_squared, ld.ci_low, ld.ci_high,
      [ld.freqs.aa, ld.freqs.ab, ld.freqs.ba, ld.freqs.bb], [x.decode('utf8') for x in ld.phase])

cdef extern from 'utils.h' namespace 'ld_estimator':
    bool is_monomorphic(vector[vector[string]]) except +

def _is_monomorphic(var):
  return is_monomorphic(to_bytes(var))
