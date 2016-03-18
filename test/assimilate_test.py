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
		now = time()
		print "running assimilate took %s seconds" % str(now-before)
		ng.delete()
		print "running delete took %s seconds"% str(time()-now)

	def test_slug_twins(self):
		# when they have very similar description fields and have trouble making unique genome names.
		ng = assimilate_from_ncbi(self.testgbfile)
		ng2 = assimilate_from_ncbi(self.testgbfile)
		self.assertNotEqual(ng.is_real(),ng2.is_real())
		ng.delete()
		ng2.delete()



if __name__ == "__main__":
	ut.main()
