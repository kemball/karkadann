from puremedian import *
from database import *
import unittest as ut
from time import time
from itertools import combinations

class medianTest(ut.TestCase):

	def test_runs(self):
		before = time()
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters where classification !='cf_putative' limit 2;")
			clusters = list(Cluster.get_many([x for (x,) in cur.fetchall()]))
		self.assertTrue(clusters[0].is_real(),clusters[1].is_real())
		pure_median_identity(clusters[0],clusters[1])
		print "1 run takes %s seconds " %(time()-before)

	def test_reflexive(self):
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters where classification !='cf_putative' limit 10;")
			clusters = list(Cluster.get_many([x for (x,) in cur.fetchall()]))
		for a,b in combinations(clusters,2):
			forward = pure_median_identity(a,b)
			reverse = pure_median_identity(b,a)
			self.assertEqual(forward,reverse)