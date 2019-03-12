
#include "utils.h"

namespace ld_estimator {

// template<typename T>
bool is_monomorphic(std::vector<std::vector<std::string> > var){
  // check that a variant has two or more alleles
  std::set<std::string> alleles;
  for (auto geno : var) {
    std::set<std::string> genoset(geno.begin(), geno.end());
    alleles.insert(genoset.begin(), genoset.end());
    if (alleles.size() > 1) {
      return false;
    }
  }
  return true;
}

bool lacks_haplotypes(Haps<int> counts, int unknown) {
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

// template<typename T, typename A>
std::vector<std::string> get_alleles(std::vector<std::vector<std::string> > var){
  // get the major and minor alleles for a variant
  //
  // The minor allele is the second most common allele in the population.
  // We've previously excluded non-polymorphic sites, so the variant has at
  // least two alleles.
  std::map<std::string, int> counts;
  for (auto geno : var) {
    for (auto allele : geno) {
      if (counts.find(allele) == counts.end()){
        counts[allele] = 0;
      }
      counts[allele] += 1;
    }
  }
  
  // create a vector of pairs, then sort the vector by the allele counts in
  // descending order
  std::vector<std::pair<std::string, int> > pairs;
  for (auto itr = counts.begin(); itr != counts.end(); ++itr){
      pairs.push_back(*itr);
  }
  std::sort(pairs.begin(), pairs.end(), [=](std::pair<std::string, int> a, std::pair<std::string, int> b)
  {
      return a.second > b.second;
  });
  
  return std::vector<std::string> {pairs[0].first, pairs[1].first};
}

}
