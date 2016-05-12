from karkadann.hmm import *
from karkadann.database import get_cursor, data_location
from karkadann.assimilate import assimilate_from_ncbi
import unittest as ut
from time import time


class ProfileTest(ut.TestCase):
	@classmethod
	def setUpClass(cls):
		testgbfile = os.path.join(data_location, 'test/testassem.gb')
		cls.ng = assimilate_from_ncbi(testgbfile)

	@classmethod
	def tearDownClass(cls):
		cls.ng.delete()

	def test_run(self):
		print "testing profile"
		assems = ProfileTest.ng.assemblies()
		for a in assems:
			contigs = a.contigs()
			for c in contigs:
				genes = list(c.genes())
				profile(genes, ["AfsA"])
				nhit = Hit(gene=genes[0],score=5,seq='LEET',hmm='notarealhmm')
				nhit.save()
				self.assertEqual(nhit.seq,'LEET')
				self.assertEqual(nhit.score,5)
				self.assertGreater(len(genes[0].hits()),0)
				nhit.delete()
				self.assertFalse(nhit.is_real())
			with get_cursor() as cur:
				# there are guaranteed to be some.
				# double checking that there's a hit from genes we just used is harder
				# but also there's no other way to put hits into the database... so
				cur.execute("select count(id) from hits;")
				for res in cur.fetchone():
					self.assertIsNotNone(res)

	def test_scan_assembly(self):
		print "testing scan_assembly"
		before = time()
		assem = next(ProfileTest.ng.assemblies())
		for contig in assem.contigs():
			for gene in contig.genes():
				for h in gene.hits():
					h.delete()
					self.assertFalse(h.is_real())
		now = time()
		print "wiped hits, took %s seconds" %(now-before)
		scan_assembly(next(ProfileTest.ng.assemblies()))
		print "scanned assembly, took %s seconds" % (time()-now)
		before = time()
		assem = next(ProfileTest.ng.assemblies())
		for contig in assem.contigs():
			for gene in contig.genes():
				for h in gene.hits():
					h.delete()
					self.assertFalse(h.is_real())
		now = time()
		print "wiped hits, took %s seconds" %(now-before)


class HmmTest(ut.TestCase):
	# TODO add more of these

	def test_list_hmms(self):
		print "testing list_hmms"
		self.assertTrue(list_hmms())
		self.assertIn("AfsA", list_hmms())
		print "there are %s hmms" % len(list_hmms())

