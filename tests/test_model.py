import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
from model import *

class ModelTests(unittest.TestCase):
    def test_design(self):
        for mean in [3,6,12,24]:
            for h in ['low','medium','high']:
                self.assertEqual(activity_counts(1000,mean,h).mean(),mean)
        b,z,d,_,active,ratios=generate(100,30,3,6,'high','broad_active',2,np.random.default_rng(3))
        np.testing.assert_array_equal(b.sum(1),d)
        self.assertEqual(active.sum(),20)
        self.assertTrue(np.all(ratios[active]==1))
        ws=projections(b)
        self.assertAlmostEqual(np.triu(ws['count'],1).sum(),np.sum(d*(d-1)/2))
        self.assertAlmostEqual(np.triu(ws['weighted'],1).sum(),np.sum(d[d>=2]/2))
        for w in ws.values():
            np.testing.assert_allclose(w,w.T)
            self.assertTrue(np.all(w.diagonal()==0))

    def test_matching_and_zeros(self):
        z=np.array([0,0,1,1,2,2])
        self.assertEqual(metrics(z,np.array([2,2,0,0,1,1]),3)['error'],0)
        self.assertAlmostEqual(metrics(z,np.array([-1,2,0,0,1,1]),3)['error'],1/6)
        w=np.zeros((6,6))
        self.assertTrue(np.all(spectral(w,3,1)==-1))
        self.assertTrue(np.all(bipartite_svd(np.zeros((5,6)),3,1)==-1))
        for w in projections(np.zeros((5,6))).values():
            self.assertTrue(np.isfinite(w).all())

    def test_uniform_theory(self):
        for d in [0,1,2,12,30]:
            q=d*(d-1)/(30*29)
            np.testing.assert_allclose(pair_moments(30,3,d,1.),[q,q],atol=1e-12)

    def test_theory_against_sampling(self):
        # Large collection of independent users makes exchangeable averages precise.
        b,z,d,_,_,ratios=generate(10000,30,3,6,'low','focused',4,np.random.default_rng(9))
        actual=edge_stats(projections(b)['count'],z)
        theory=theoretical_weights(30,3,d,ratios)['count']
        np.testing.assert_allclose([actual['within'],actual['between']],theory,rtol=.02)

    def test_partition_recovery(self):
        b=np.kron(np.eye(3),np.ones((10,5)))
        z=np.repeat(np.arange(3),5)
        for w in projections(b).values():
            self.assertEqual(metrics(z,spectral(w,3,5),3)['error'],0)
        self.assertEqual(metrics(z,bipartite_svd(b,3,5),3)['error'],0)

if __name__=='__main__':
    unittest.main()
