from karkadann.domains import *
from karkadann.domains import _call_mcl

from karkadann.database import *

import unittest as ut

class mclTest(ut.TestCase):
	def test_call_mcl(self):
		print "test_call_mcl"
		with get_cursor() as cur:
			cur.execute("select id from genes limit 1000;")
			records = [Gene(db_id=x).record for (x,) in cur.fetchall()]
			_call_mcl(records)

	def test_locking(self):
		print"test_locking"
		import multiprocessing as mp
		with get_cursor() as cur:
			cur.execute("select id from genes limit 1000;")
			records = [Gene(db_id=x).record for (x,) in cur.fetchall()]
			p1 = mp.Process(target=_call_mcl,args=(records,))
			cur.execute("select id from genes limit 1000;")
			records2 = [Gene(db_id=x).record for (x,) in cur.fetchall()]
			p2 = mp.Process(target=_call_mcl,args=(records2,))
		p1.start()
		p2.start()
		p1.join()
		p2.join()
		# how to check if either broke? I guess if there's database sharing at least one won't complete
		# mcl is fragile in that way.

	def test_parsing(self):
		print "test_parsing"
		with get_cursor() as cur:
			cur.execute("select id from genes limit 1000;")
			records = [Gene(db_id=x).record for (x,) in cur.fetchall()]
			gtext = _call_mcl(records)
			parse_groups(gtext)
			# i guess it doesn't crash at least.

	def test_start_batch(self):
		print "test_batch"
		batchid = start_batch()
		testid = most_recent_batch()
		self.assertEqual(batchid,testid)

	def test_holistic(self):
		print "test all the things"
		with get_cursor() as cur:
			cur.execute("select id from contigs order by length(sequence) desc limit 1 ;")
			andre = Contig(db_id=cur.fetchone()[0])
			print "Andre secured"
		nb = start_batch()
		try:
			assign_groups(list(andre.genes())[:2000])
		finally:
			with get_cursor() as cur:
				cur.execute("delete from orthomcl_batches where id = %s",(nb,))