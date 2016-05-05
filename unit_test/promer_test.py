# TODO write some tests for promer

from karkadann.promer import *
from karkadann.database import *
from karkadann.assimilate import assimilate_from_ncbi
from karkadann.hmm import scan_assembly
from karkadann.cluster_call import call_clusters
import unittest as ut


class PromerTest(ut.TestCase):
	#TODO a more complete test case for this is necessary

	@classmethod
	def setUpClass(cls):
		with get_cursor() as cur:
			cur.execute("select count(*) from clusters;")
			if not cur.fetchone()[0]:
				testgbfile = os.path.join(data_location, 'test/testassem.gb')
				cls.ng = assimilate_from_ncbi(testgbfile)
				cls.assem = cls.ng.assemblies().next()
				scan_assembly(cls.assem)
				for contig in cls.assem.contigs():
					call_clusters(contig)

	def test_identity(self):
		with get_cursor() as cur:
			cur.execute("select id from clusters limit 1;")
			cluster = Cluster.get(cur.fetchone()[0])
		score = promer_score(cluster,cluster)
		self.assertEqual(score,1)


	def test_promer_runs(self):
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters limit 2;")
			one_id = cur.fetchone()[0]
			two_id = cur.fetchone()[0]
		one = Cluster.get(db_id=one_id)
		two = Cluster.get(db_id=two_id)
		score = promer_score(one, two)
		self.assertGreaterEqual(score, 0)
		self.assertLessEqual(score, 1)
		score = promer_score(one, two)
		self.assertGreaterEqual(score, 0)
		self.assertLessEqual(score, 1)
