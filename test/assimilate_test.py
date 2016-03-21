from karkadann.assimilate import *
from karkadann.assimilate import _slug
from karkadann.database import data_location
from Bio import SeqIO
import os

import unittest as ut


class AssimilateTest(ut.TestCase):
	testgbfile = os.path.join(data_location, 'test/testassem.gb')

	def test_slug_basic(self):
		recs = list(SeqIO.parse(self.testgbfile, 'genbank'))
		desc = recs[0].description
		self.assertNotEqual(_slug(desc, 3), _slug(desc, 4))

	def test_standardize(self):
		recs = list(SeqIO.parse(self.testgbfile, 'genbank'))
		self.assertNotIn('translation', recs[0].features[4].qualifiers.keys())
		make_standard(recs)
		self.assertIn('translation', recs[0].features[4].qualifiers.keys())

	def test_assimilate(self):
		from time import time
		before = time()
		ng = assimilate_from_ncbi(self.testgbfile)
		assem = next(ng.assemblies())
		from karkadann.database import get_cursor

		now = time()
		print "running assimilate took %s seconds" % str(now-before)
		with get_cursor() as cur:
			cur.execute("select accession from assemblies where id = %s;",(assem.is_real()))
			self.assertEqual(cur.fetchone()[0],"GCF_000576425.1")
		with get_cursor() as cur:
			cur.execute("select accession from contigs where accession=%s;",("NZ_KK070022.1",))
			self.assertIsNotNone(cur.fetchone())
		with get_cursor() as cur:
			cur.execute("select accession from genes where accession=%s;",("WP_038374786.1",))
			self.assertIsNotNone(cur.fetchone())
		ng.delete()
		print "running delete took %s seconds" % str(time()-now)

	def test_slug_twins(self):
		# when they have very similar description fields and have trouble making unique genome names.
		ng = assimilate_from_ncbi(self.testgbfile)
		ng2 = assimilate_from_ncbi(self.testgbfile)
		self.assertNotEqual(ng.is_real(),ng2.is_real())
		self.assertNotEqual(ng._name,ng2._name)
		ng.delete()
		ng2.delete()


if __name__ == "__main__":
	ut.main()
