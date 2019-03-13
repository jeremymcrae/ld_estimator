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

cdef extern from 'ld.h' namespace 'ld_estimator':
    Linkage pairwise(vector[vector[string]], vector[vector[string]], vector[bool]) except +

def to_bytes(var):
    return [list(map(str.encode, x)) for x in var]

def pairwise_ld(var1, var2, vector[bool] ploidy):
    try:
        ld = pairwise(to_bytes(var1), to_bytes(var2), ploidy)
    except ValueError:
        return None
    
    return LD(ld.dprime, ld.loglikelihood, ld.r_squared, ld.ci_low, ld.ci_high,
      [ld.freqs.aa, ld.freqs.ab, ld.freqs.ba, ld.freqs.bb])
