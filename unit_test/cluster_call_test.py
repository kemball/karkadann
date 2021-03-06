from karkadann.cluster_call import _classify
from karkadann.cluster_call import *
from karkadann.database import get_cursor
import unittest as ut
from time import time


class ClassifyTest(ut.TestCase):
	@classmethod
	def setUpClass(cls):
		from karkadann.assimilate import assimilate_from_ncbi
		from karkadann.hmm import scan_assembly
		from karkadann.database import data_location
		print "assimilating for classification testing"
		before = time()
		cls.ng = assimilate_from_ncbi(data_location+'/'+'test/test2.gb')
		cls.assem = next(cls.ng.assemblies())
		cls.contigs = list(cls.assem.contigs())
		print "assimilation took %d seconds" % (time()-before)
		before = time()
		scan_assembly(cls.assem)
		print "scanning an assembly takes %s seconds" % (time()-before)

	@classmethod
	def tearDownClass(cls):
		cls.ng.delete()

	def test_classify(self):
		gene = list(ClassifyTest.contigs[-1].genes())[-1]
		before = time()
		_classify(gene)
		print "classifying a single gene takes %s seconds"%(time()-before)

	def test_call_clusters(self):
		for contig in ClassifyTest.contigs:
			before = time()
			call_clusters(contig)
			print "calling clusters in a single contig takes %s seconds " % (time() - before)
		for contig in ClassifyTest.contigs:
			with get_cursor() as cur:
				cur.execute("select distinct clusters.id as clust_id,(select count(*) from clusters where id=clust_id) from clusters join genes on genes.id=clusters.gene join contigs on genes.contig=contigs.id where contigs.id=%s;",(contig.is_real(),))
				for (clust_id,sum) in cur.fetchall():
					self.assertGreater(sum,1)
			print "good, no singleton clusters."






