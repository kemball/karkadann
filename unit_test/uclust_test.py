import unittest as ut
import karkadann
from karkadann.uclust import *
from karkadann.database import get_cursor,Gene,domain_max


class uclustTest(ut.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def test_call_uclust(self):
		print "testing uclust_all"
		with get_cursor() as cur:
			cur.execute("select distinct id from genes limit 500;")
			results = cur.fetchall()
			genes = Gene.get_many([x for (x,) in results])
		hitlist = list(itertools.chain(*[g.hits() for g in genes]))
		if not len(hitlist):
			return
		d = uclust_all(hitlist,.5)
		self.assertEqual(len(hitlist),sum([len(v) for v in d.values()]))
		self.assertGreater(len(d.keys()),0)
		for v in d.values():
			self.assertTrue(v)
			for hitid in v:
				self.assertTrue(hitid.is_real())

	def test_num_clusters(self):
		print "testing that no hits are lost"
		with get_cursor() as cur:
			cur.execute("select distinct id from genes limit 1000;")
			results = cur.fetchall()
			genes = Gene.get_many([x for (x,) in results])
		hitlist = list(itertools.chain(*[g.hits() for g in genes]))
		d = uclust_all(hitlist,.5)
		self.assertEqual(len(hitlist),sum([len(v) for v in d.values()]))


	def test_set_of_domains(self):
		print "testing domain set"
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters limit 2;")
			clusters = [Cluster.get(db_id) for (db_id,) in cur.fetchall()]
		if not clusters:
			return
		# TODO make this test something.

	def test_doroghazi_uclust(self):
		print "testing domain_max calculation"
		calc_domain_max("nrps")
		calc_domain_max("PKS_I")
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters where classification='nrps' limit 2;")
			clusters = [Cluster.get(db_id) for (db_id,) in cur.fetchall()]
			if not clusters:
				return
			self.assertGreaterEqual(domain_max(clusters[0],clusters[1]),0)
			self.assertEqual(domain_max(*clusters),domain_max(clusters[1],clusters[0]))
			self.assertEqual(domain_max(clusters[0],clusters[0]),1.0)




