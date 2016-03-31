# TODO write some tests for promer

from karkadann.promer import *
from karkadann.database import *
import unittest as ut


class PromerTest(ut.TestCase):
	pass

	def test_promer_runs(self):
		one = Cluster(db_id="nrpsq1")
		two = Cluster(db_id="nrpsq2")
		score = promer_score(one, two)
		self.assertGreaterEqual(score, 0)
		self.assertLessEqual(score, 1)
