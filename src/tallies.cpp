
#include "tallies.h"

namespace ld_estimator {

// template<typename T>
std::pair<Haps<int>, int> tally_haplotypes(std::vector<std::vector<std::string> > var1,
    std::vector<std::vector<std::string> > var2, std::vector<bool> ploidy) {
  // tally known haplotypes of known and unknown phase
  //
  // Args:
  //     var1: list of genotypes for first variant
  //     var2: list of genotypes for second variant, ordered as per var1
  //     ploidy: list of ploidy states, ordered as per var1
  int unknown = 0;
  std::vector<std::string> alleles = get_alleles(var1);
  std::string major_1 = alleles[0];
  std::string minor_1 = alleles[1];
  alleles = get_alleles(var2);
  std::string major_2 = alleles[0];
  std::string minor_2 = alleles[1];
  
  std::map<std::string, std::map<std::string, int>> slots;
  slots[major_1][major_2] = 0;
  slots[major_1][minor_2] = 0;
  slots[minor_1][major_2] = 0;
  slots[minor_1][minor_2] = 0;

  // iterate through all chromosomes in dataset
  int size = var1.size();
  for (int i=0; i < size; i++) {
    std::string a1 = var1[i][0];
    std::string a2 = var1[i][1];
    std::string b1 = var2[i][0];
    std::string b2 = var2[i][1];
    bool is_haploid = ploidy[i];
    if (!is_haploid) {
      if (a1 == "None" or a2 == "None" or b1 == "None" or b2 == "None") {
        // skip missing data
        continue;
      }
      else if (a1 != a2 and b1 != b2) {
        // both genotypes are heterozygous
        unknown += 1;
      }
      else if (a1 != a2) {
        // b alleles are homozygous
        slots[major_1][b1] += 1;
        slots[minor_1][b1] += 1;
      }
      else if (b1 != b2) {
        // a alleles are homozygous
        slots[a1][major_2] += 1;
        slots[a1][minor_2] += 1;
      }
      else {
        // both are homozygous
        slots[a1][b1] += 2;
      }
    }
    else {
      // include haploid chromosomes
      if (a1 != "None" and a2 != "None") {
        slots[a1][b1] += 1;
      }
    }
  }

  int aa = slots[major_1][major_2];
  int ab = slots[major_1][minor_2];
  int ba = slots[minor_1][major_2];
  int bb = slots[minor_1][minor_2];
  Haps<int> counts = Haps<int> {aa, ab, ba, bb};

  return std::pair<Haps<int>, int> (counts, unknown);
}

}
