from karkadann.cluster_call import _classify
from karkadann.cluster_call import *
import unittest as ut
from time import time


class ClassifyTest(ut.TestCase):
	@classmethod
	def setUpClass(cls):
		from assimilate import assimilate_from_ncbi
		from hmm import scan_assembly
		print "assimilating for classification testing"
		before = time()
		cls.ng = assimilate_from_ncbi('../test/testassem.gb')
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
		before = time()
		call_clusters(ClassifyTest.contigs[-1])
		print "calling clusters in a single contig takes %s seconds "% (time()-before)


