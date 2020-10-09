# cython: language_level=3, boundscheck=False

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
    Linkage pairwise(vector[string] & var_a1, vector[string] & var_a2,
        vector[string] & var_b1, vector[string] & var_b2, vector[bool] &) except +

cdef to_bytes(var):
    return [[str(a).encode('utf8') if a is not None else b'' for a in g] for g in var]

def pairwise_ld(var1, var2, vector[bool] ploidy):
    # only convert to allele vectors to bytes if not done already
    if not isinstance(var1[0][0], bytes):
        var1 = to_bytes(var1)
    if not isinstance(var2[0][0], bytes):
        var2 = to_bytes(var2)
    
    if len(var1) != len(var2):
      raise ValueError('genotypes lists need to be the same length')
    if len(var1) == 0:
      raise ValueError('zero genotypes in lists supplied for LD')

    cdef vector[string] var_a1
    cdef vector[string] var_a2
    cdef vector[string] var_b1
    cdef vector[string] var_b2
    
    var_a1.resize(len(var1))
    var_a2.resize(len(var1))
    var_b1.resize(len(var1))
    var_b2.resize(len(var1))
    
    for i, x in enumerate(var1):
        var_a1[i] = x[0] if x[0] is not None else b''
        var_a2[i] = x[-1] if x[-1] is not None else b''

    for i, x in enumerate(var2):
        var_b1[i] = x[0] if x[0] is not None else b''
        var_b2[i] = x[-1] if x[-1] is not None else b''
    
    try:
        ld = pairwise(var_a1, var_a2, var_b1, var_b2, ploidy)
    except ValueError:
        return None

    return LD(ld.dprime, ld.loglikelihood, ld.r_squared,
      [ld.freqs.aa, ld.freqs.ab, ld.freqs.ba, ld.freqs.bb], [ld.phase.var1_allele.decode('utf8'), ld.phase.var2_allele.decode('utf8')])

cdef extern from 'utils.h' namespace 'ld_estimator':
    bool is_monomorphic(vector[string] &, vector[string] &) except +

def _is_monomorphic(var):
    cdef vector[string] * a1 = new vector[string](len(var))
    cdef vector[string] * a2 = new vector[string](len(var))
    for i, x in enumerate(var):
        a1[i] = <string> str(x[0]).encode('utf8') if x[0] is not None else b''
        a2[i] = <string> str(x[-1]).encode('utf8') if x[-1] is not None else b''
    
    return is_monomorphic(deref(a1), deref(a2))
