import unittest as ut
import karkadann
from karkadann.uclust import *
from karkadann.database import get_cursor,Gene


class uclustTest(ut.TestCase):

	@classmethod
	def setUpClass(cls):
		pass

	@classmethod
	def tearDownClass(cls):
		pass

	def test_call_uclust(self):
		with get_cursor() as cur:
			cur.execute("select distinct id from genes limit 500;")
			results = cur.fetchall()
			genes = Gene.get_many([x for (x,) in results])
		d = uclust_all(genes,.5)
		self.assertGreater(len(d.keys()),0)
		for v in d.values():
			self.assertTrue(v)
			for hitid in v:
				self.assertTrue(Hit.get(hitid).is_real())

	def test_parse_ucclust(self):
		pass


