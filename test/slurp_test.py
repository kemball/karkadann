from context import karkadann
from karkadann.slurp import *
from karkadann.database import *
import os
from Bio import SeqIO
import unittest as ut 

class GenbankSlurp(ut.TestCase):
	def test_annotate(self):
		records = SeqIO.parse(os.path.join(data_location,'test/testassem.gb'),'genbank')
		records = list(records)
		newassem = annotate_save(os.path.join(data_location,'test/testassem.gb'),'testingslurp')
		with get_cursor() as cur:
			gid = get_genome('testingslurp')
			cur.execute("select id from assemblies where genome_id = %s;",(gid,))
			assid = cur.fetchone()[0]
			self.assertEqual(assid,newassem)
			cur.execute("select count(*) from contigs where assembly_id=%s;",(assid,))
			num_contigs = cur.fetchone()[0]
		self.assertEqual(num_contigs,len(records))

		with get_cursor() as cur:
			gid = get_genome('testingslurp')
			cur.execute("select id from assemblies where genome_id = %s;",(gid,))
			assemid = cur.fetchone()[0]
			self.assertEqual(assemid,newassem)
			contig_ids = get_contigs(assemid)
			cur.executemany("delete from genes where contig = %s;",[(x) for x in contig_ids])
			cur.executemany("delete from contigs where id = %s;",[(x) for x in get_contigs(assemid)])
			cur.execute("delete from assemblies where id = %s",(assemid,))
			cur.execute("delete from genus_species where genome_id=%s;",(gid,))
			cur.execute("delete from genomes where id = %s;",(gid,))








if __name__ == '__main__':
	ut.main()