''' unit tests for LD class
'''

import unittest

from ld_estimator import pairwise_ld

class TestLDEstimator(unittest.TestCase):
    ''' class to unit test LD calculations
    '''
    def test_pairwise_ld(self):
        ''' check pairwise_ld works correctly
        '''
        var1 = [('A', 'A'), ('A', 'A'), ('A', 'G'), ('G', 'A'), ('G', 'G'),
                ('A', 'G'), ('A', 'A'), ('A', 'A'), ('G', 'G')]
        var2 = [('A', 'A'), ('A', 'A'), ('A', 'G'), ('G', 'A'), ('G', 'G'),
                ('A', 'G'), ('A', 'A'), ('G', 'G'), ('G', 'G')]
        is_haploid = [False, False, False, False, False, False, False, True, True, True]

        linkage = pairwise_ld(var1, var2, is_haploid)
        self.assertEqual(linkage.dprime, 0.9999999995259259)
        self.assertEqual(linkage.r_squared, 0.7714285713447618)

    def test_pairwise_ld_low_dprime(self):
        ''' check pairwise_ld works correctly for variants with low D'
        '''
        var1 = [('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'),
                ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'),
                ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'),
                ('G', 'G'), ('G', 'A'), ('G', 'G'), ('G', 'A'), ('G', 'G'),
                ('G', 'G'), ('G', 'G'), ('G', 'G')]

        var2 = [('G', 'G'), ('G', 'G'), ('A', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'),
                ('G', 'A'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('A', 'G'), ('G', 'G'),
                ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('A', 'G'), ('G', 'G'),
                ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'A')]

        is_haploid = [False] * len(var1)
        linkage = pairwise_ld(var1, var2, is_haploid)
        self.assertEqual(linkage.dprime, 0.37448437109790317)
        self.assertEqual(linkage.r_squared, 0.05227073010963888)

    def test_pairwise_ld_monomorphic(self):
        ''' test we don't compute LD with monomorphic markers
        '''
        var1 = [('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G')]
        var2 = [('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G')]
        is_haploid = [False] * len(var1)
        self.assertIsNone(pairwise_ld(var1, var2, is_haploid))

        # check it gives none even if only one marker is monomorphic
        var3 = [('G', 'A'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G'), ('G', 'G')]
        self.assertIsNone(pairwise_ld(var1, var3, is_haploid))

    def test_lacks_haplotypes(self):
        ''' test pairwise_ld when missing some haplotypes
        '''
        var1 = [('0', '0'), ('0', '0'), ('0', '0'), ('1', '1'), (None, None)]
        var2 = [('0', '0'), ('0', '0'), ('0', '0'), (None, None), ('1', '1')]
        is_haploid = [False] * len(var1)
        linkage = pairwise_ld(var1, var2, is_haploid)
        self.assertEqual(linkage.dprime, 1)
        self.assertEqual(linkage.r_squared, 0)
        self.assertEqual(linkage.freqs, [0, 0, 0, 0])
