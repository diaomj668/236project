import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import unittest
import numpy as np
from model import generate_separation, append_background, projections, spectral, metrics
from supplement import seed, DESIGN
from threadpoolctl import threadpool_limits

class SupplementTests(unittest.TestCase):
    def setUp(self):
        self.signal,self.z,self.order=generate_separation(100,30,3,6,3,20,24,np.random.default_rng(1),np.random.default_rng(2))
    def test_fixed_signal_and_nested_background(self):
        b=append_background(self.signal,self.order,10,12)
        larger=append_background(self.signal,self.order,20,24)
        np.testing.assert_array_equal(b[:100],self.signal)
        np.testing.assert_array_equal(larger[:100],self.signal)
        np.testing.assert_array_equal(b[:100].sum(1),6)
        np.testing.assert_array_equal(b[100:].sum(1),12)
        self.assertTrue(np.all(b[100:]<=larger[100:110]))
        self.assertEqual(append_background(self.signal,self.order,0,6).shape,(100,30))
    def test_constant_activity_equivalence(self):
        b=append_background(self.signal,self.order,20,6)
        w=projections(b)
        np.testing.assert_allclose(w['weighted'],w['count']/5)
        def norm(a):
            d=a.sum(1);return a/np.sqrt(d[:,None]*d[None,:])
        np.testing.assert_allclose(norm(w['count']),norm(w['weighted']),atol=1e-12)
        with threadpool_limits(limits=1):
            c=spectral(w['count'],3,23);v=spectral(w['weighted'],3,23)
        self.assertEqual(metrics(self.z,c,3)['error'],metrics(self.z,v,3)['error'])
    def test_seed_namespaces(self):
        self.assertEqual(seed(23604,1,2),seed(23604,1,2))
        self.assertNotEqual(seed(23604,1,2),seed(23605,1,2))
        self.assertEqual(DESIGN['separation']['noise_users'],[0,50,200])

if __name__=='__main__':unittest.main()
