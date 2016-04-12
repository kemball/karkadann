from karkadann.domains import *
from karkadann.domains import _call_mcl

from karkadann.database import *

import unittest as ut
from itertools import chain

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

	def test_save_orthogroup(self):
		print 'test_save_orthogroup'
		with get_cursor() as cur:
			cur.execute("select id from genes limit 1;")
			gene = Gene(db_id=cur.fetchone()[0])
			save_orthogroup(gene,'ortho1337')
		self.assertEqual(gene.orthogroup,"ortho1337")


	def test_holistic(self):
		print "test all the things"
		with get_cursor() as cur:
			cur.execute("select id from contigs order by length(sequence) desc limit 1 ;")
			andre = Contig(db_id=cur.fetchone()[0])
			print "Andre secured"
		nb = start_batch()
		try:
			assign_groups(list(andre.genes())[:2000])
			print "done assigning groups"
		finally:
			with get_cursor() as cur:
				cur.execute("delete from orthomcl_batches where id = %s",(nb,))

	def test_score(self):
		print "testing domain scoring"
		with get_cursor() as cur:
			# PKS are very common but come in three kinds.
			cur.execute("select  id from clusters limit 5 where classification like 'PKS%';")
			clusters = [Cluster(db_id = x) for (x,) in cur.fetchall()]
		nb = start_batch()
		try:
			genes = chain(*[c.gene_list() for c in clusters])
			assign_groups(genes)
			clusta = clusters[0]
			clustb = clusters[1]
			self.assertEqual(domain_score(clusta,clustb),domain_score(clustb,clusta))
			self.assertGreaterEqual(domain_score(clusters[1],clusters[2]),0)
			self.assertLessEqual(domain_score(clusters[3],clusters[4]),1)
		finally:
			with get_cursor() as cur:
				cur.execute("delete from orthomcl_batches where id = %s",(nb,))