
#include <iostream>
#include <array>
#include <cstring>

#include "tallies.h"

namespace ld_estimator {

// template<typename T>
std::pair<Haps<int>, int> tally_haplotypes(char ** a1s, char ** a2s,
    char ** b1s, char ** b2s, int size, std::vector<bool> & ploidy,
    const char * major_1, const char * minor_1, const char * major_2, const char * minor_2) {
  // tally known haplotypes of known and unknown phase
  //
  // Args:
  //     var1: list of genotypes for first variant
  //     var2: list of genotypes for second variant, ordered as per var1
  //     ploidy: list of ploidy states, ordered as per var1
  int unknown = 0;
  std::map<const char *, std::map<const char *, int>> slots;
  slots[major_1][major_2] = 0;
  slots[major_1][minor_2] = 0;
  slots[minor_1][major_2] = 0;
  slots[minor_1][minor_2] = 0;
  
  char * a1;
  char * a2;
  char * b1;
  char * b2;
  // bool is_haploid;
  // iterate through all chromosomes in dataset
  for (int i=0; i < size; i++) {
    a1 = a1s[i];
    b1 = b1s[i];
    if (!ploidy[i]) {
      a2 = a2s[i];
      b2 = b2s[i];
      if (std::strlen(a1) == 0 || std::strlen(a2) == 0 || std::strlen(b1) == 0 || std::strlen(b2) == 0) {
        // skip missing data
        continue;
      }
      else if ((std::strcmp(a1, a2) != 0) & (std::strcmp(b1, b2) != 0)) {
        // both genotypes are heterozygous
        unknown += 1;
      }
      else if (std::strcmp(a1, a2) != 0) {
        // b alleles are homozygous
        slots[major_1][b1] += 1;
        slots[minor_1][b1] += 1;
      }
      else if (std::strcmp(b1, b2) != 0) {
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
      if ((std::strlen(a1) != 0) & (std::strlen(a2) != 0) ) {
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
