
import unittest

from ld_estimator import estimate_ld

class TestLDEstimator(unittest.TestCase):
    
    def test_estimate_ld(self):
        ''' check estimate_ld works correctly
        '''
        var1 = [(0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (0, 1), (0, 0), (0, 0), (1, 1)]
        var2 = [(0, 0), (0, 0), (0, 1), (1, 0), (1, 1), (0, 1), (0, 0), (1, 1), (1, 1)]
        is_haploid = [False, False, False, False, False, False, False, True, True, True]
        
        ld = estimate_ld(var1, var2, is_haploid)
        self.assertEqual(ld.dprime, 0.9999999995259259)
        self.assertEqual(ld.r_squared, 0.7714285713447618)
    
    def test_estimate_ld_low_dprime(self):
        ''' check estimate_ld works correctly for variants with low D'
        '''
        var1 = [('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G'),
            ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G'),
            ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','A'), ('G','G'),
            ('G','A'), ('G','G'), ('G','G'), ('G','G'), ('G','G')]

        var2 = [('G','G'), ('G','G'), ('A','G'), ('G','G'), ('G','G'), ('G','G'),
            ('G','A'), ('G','G'), ('G','G'), ('G','G'), ('A','G'), ('G','G'),
            ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('A','G'), ('G','G'),
            ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','A')]

        is_haploid = [False] * len(var1)
        ld = estimate_ld(var1, var2, is_haploid)
        self.assertEqual(ld.dprime, 0.37448437109790317)
        self.assertEqual(ld.r_squared, 0.05227073010963888)
    
    def test_estimate_ld_monomorphic(self):
        ''' test we don't compute LD with monomorphic markers
        '''
        var1 = [('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G')]
        var2 = [('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G')]
        is_haploid = [False] * len(var1)
        self.assertIsNone(estimate_ld(var1, var2, is_haploid))

        # check it gives none even if only one marker is monomorphic
        var3 = [('G','A'), ('G','G'), ('G','G'), ('G','G'), ('G','G'), ('G','G')]
        self.assertIsNone(estimate_ld(var1, var3, is_haploid))
