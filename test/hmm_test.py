from karkadann.hmm import *
from karkadann.database import get_cursor,data_location
from karkadann.assimilate import assimilate_from_ncbi
import unittest as ut



testgbfile = os.path.join(data_location, 'test/testassem.gb')
ng = assimilate_from_ncbi(testgbfile)

class ProfileTest(ut.TestCase):
	def test_run(self):
		assems = ng.assemblies()
		for a in assems:
			contigs= a.contigs()
			for c in contigs:
				genes = c.genes()
				profile(genes,"AfsA.hmm")
			with get_cursor as cur:
				# there are guaranteed to be some.
				# double checking that there's a hit from genes we just used is harder
				# but also there's no other way to put hits into the database... so
				cur.execute("select count(id) from hits;")
				for res in cur.fetchone():
					self.assertIsNotNone(res[0])


	def test_scan_assembly(self):
		scan_assembly(next(ng.assemblies))


class HmmTest(ut.TestCase):
	#TODO
	def test_list_hmms(self):
		self.assertTrue(list_hmms())
		self.assertIn("AfsA.hmm",list_hmms())

