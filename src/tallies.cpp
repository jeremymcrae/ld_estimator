
#include <iostream>
#include <array>

#include "tallies.h"

namespace ld_estimator {

// template<typename T>
std::pair<Haps<int>, int> tally_haplotypes(std::vector<std::string> & a1s, std::vector<std::string> & a2s,
    std::vector<std::string> & b1s, std::vector<std::string> & b2s, std::vector<bool> & ploidy,
    std::string major_1, std::string minor_1, std::string major_2, std::string minor_2) {
  // tally known haplotypes of known and unknown phase
  //
  // Args:
  //     var1: list of genotypes for first variant
  //     var2: list of genotypes for second variant, ordered as per var1
  //     ploidy: list of ploidy states, ordered as per var1
  int unknown = 0;
  std::map<std::string, std::map<std::string, int>> slots;
  slots[major_1][major_2] = 0;
  slots[major_1][minor_2] = 0;
  slots[minor_1][major_2] = 0;
  slots[minor_1][minor_2] = 0;
  
  std::string a1;
  std::string a2;
  std::string b1;
  std::string b2;
  // bool is_haploid;
  // iterate through all chromosomes in dataset
  int size = a1s.size();
  for (int i=0; i < size; i++) {
    a1 = a1s[i];
    b1 = b1s[i];
    if (!ploidy[i]) {
      a2 = a2s[i];
      b2 = b2s[i];
      if (a1.empty() || a2.empty() || b1.empty() || b2.empty()) {
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
      a2 = a2s[i];
      // include haploid chromosomes
      if (a1.empty() & a2.empty()) {
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
