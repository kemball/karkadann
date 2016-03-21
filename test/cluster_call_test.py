from karkadann.cluster_call import _classify
from karkadann.cluster_call import *
import unittest as ut


class classifyTest(ut.TestCase):
	@classmethod
	def setUpClass(cls):
		from assimilate import assimilate_from_ncbi
		from hmm import scan_assembly
		cls.ng = assimilate_from_ncbi('test/testassem.gb')
		scan_assembly(next(ng.assemblies()))

	@classmethod
	def tearDownClass(cls):
		cls.ng.delete()

	def test_classify(self):
		gene = next(next(next(ng.assemblies()).contigs()).genes())
		_classify(gene)




class clusterTest(ut.TestCase):
	@classmethod
	def setUpClass(cls):
		from assimilate import assimilate_from_ncbi
		from hmm import scan_assembly
		cls.ng = assimilate_from_ncbi('test/testassem.gb')
		cls.assem = next(ng.assemblies())
		cls.contigs = list(cls.assem.contigs())
		scan_assembly(cls.assem)

	@classmethod
	def tearDownClass(cls):
		cls.ng.delete()

	def test_call_clusters(self):
		call_clusters(cls.contigs[-1])