
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
