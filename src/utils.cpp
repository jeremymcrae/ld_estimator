
#include <set>

#include "utils.h"

namespace ld_estimator {

// template<typename T>
bool is_monomorphic(char ** a1, char ** a2, int size){
  // check that a variant has two or more alleles
  std::set<const char *> alleles;
  for (int i=0; i < size; i++) {
    alleles.insert(a1[i]);
    alleles.insert(a2[i]);
    if (alleles.size() > 1) {
      return false;
    }
  }
  return true;
}

bool lacks_haplotypes(Haps<int> & counts, int unknown) {
  /* check if the two variants lack some haplotypes

  Even though monomorphic markers have been excluded, we need to check that
  different haplotypes are present. It is possible to lack certain haplotypes
  if the genotypes with the minor allele are paired to genotypes with
  genotypes with missing data.

  e.g.
  from ld_estimator.tallies import count_haplotypes
  var1 = [(0,0), (0,0), (0,0), (1,1), (None,None)]
  var2 = [(0,0), (0,0), (0,0), (None,None), (1,1)]
  ploidy = [False] * len(var1)
  counts, unknown = count_haplotypes(var1, var2, ploidy)
  lacks_haplotypes(counts, unknown)
  */
  
  int r1 = counts.aa + counts.ab;
  int r2 = counts.ba + counts.bb;
  int c1 = counts.aa + counts.ba;
  int c2 = counts.ab + counts.bb;
  return (r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0) and unknown == 0;
}

// Comparison function for sorting the set by decreasing order of its pair's
// second value
struct revcomp {
  template <typename T>

  // Comparator function
  bool operator()(const T &l, const T &r) const
  {
    if (l.second != r.second)  {
      return l.second > r.second;
    }
    return l.first < r.first;
  }
};

// template<typename T, typename A>
void get_alleles(char ** a1, char ** a2, int size, std::vector<bool> & haploid, std::array<const char *, 2> & alleles){
  // get the major and minor alleles for a variant
  //
  // The minor allele is the second most common allele in the population.
  // We've previously excluded non-polymorphic sites, so the variant has at
  // least two alleles.
  std::map<const char *, int> counts;
  for (int i=0; i < size; i++) {
      counts[a1[i]] += 1;
  }
  for (int i=0; i < size; i++) {
      if (haploid[i]) {
        // if the sample is haploid, don't consider the second allele
        continue;
      }
      counts[a2[i]] += 1;
  }

  std::set<std::pair<const char *, int>, revcomp> pairs(counts.begin(), counts.end());

  auto it = pairs.begin();
  alleles[0] = (*it).first;
  std::advance(it, 1);
  alleles[1] = (*it).first;
}

}
