from context import karkadann
from karkadann.database import *
import unittest as ut 

class FileTest(ut.TestCase):
	import os
	def test_data_location(self):
		self.assertTrue(os.access(data_location,os.R_OK))
	def test_gb_location(self):
		self.assertTrue(os.access(gb_location,os.R_OK))
	def test_hmm_location(self):
		self.assertTrue(os.access(hmm_location,os.R_OK))




class GenomeTest(ut.TestCase):
	def test_make_genome(self):
		with get_cursor() as cursorone:
			with get_cursor() as cursortwo:
				self.assertIsNot( cursorone , cursortwo)

		testid = make_genome("thisisnotarealname")
		self.assertTrue(testid)
		with get_cursor() as cur:
			cur.execute("delete from genomes where id = %s;",(testid,))

	def test_genome_db(self):
		testid = make_genome("thisisnotarealname")
		retrievedid = get_genome("thisisnotarealname")
		self.assertEqual(testid,retrievedid)
		with get_cursor() as cur:
			cur.execute("delete from genomes where id = %s",(testid,))
			#so that shouldn't commit until it falls otu of scope, right?
			self.assertTrue( get_genome("thisisnotarealname"))
		self.assertFalse(get_genome("thisisnotarealname"))

	def test_implementation(self):
		with get_cursor() as cursorone:
			cursorone.execute("insert into genomes (name,id) values (%s,%s);",("thisisnotarealname",42))
			with get_cursor() as cursortwo:
				self.assertFalse( get_genome("thisisnotarealname")) 
		self.assertTrue(get_genome("thisisnotarealname"))
		with get_cursor() as cursor:
			cursor.execute("delete from genomes where id=%s",(get_genome("thisisnotarealname"),))



class AssemblyTest(ut.TestCase):
	#this screams for setup and teardown methods
	#and yet here we are
	records = SeqIO.parse(os.path.join(data_location,'test/testassem.gb'),'genbank')

	def test_make_assembly(self):
		records = list(self.records)
		gid = make_genome('test')
		aid = make_assembly(records,gid)
		with get_cursor() as cur:
			cur.execute('delete from assemblies where id=%s;',(aid,))
			cur.execute('delete from genomes where id =%s;',(gid,))

	def test_gb_record(self):
		records = list(self.records)
		gid = make_genome("foobarbazqux")
		newassem = make_assembly(records,gid)
		with get_cursor() as cur:
			cur.execute("select gb_record from assemblies where id = %s;",(newassem,))
			name = cur.fetchone()[0]
			self.assertTrue(os.path.exists(os.path.join(gb_location,name)))
		with get_cursor() as cur:
			cur.execute("select genomes.id,genomes.name from genomes\
														join assemblies \
															on assemblies.genome_id=genomes.id\
														where assemblies.id=%s;",(newassem,))
			fid,fname = cur.fetchone()
			self.assertEqual(gid,fid)
			self.assertEqual("foobarbazqux",fname)
		roundaboutrecord = list(read_record(newassem))
		assert(len(roundaboutrecord)==len(records))
		assert(roundaboutrecord[0].id == records[0].id)
		os.remove(os.path.join(gb_location,name))
		with get_cursor() as cur:
			cur.execute('delete from assemblies where id=%s;',(assemid,))
			cur.execute('delete from genomes where id =%s;',(gid,))





class HmmTest(ut.TestCase):
	def test_import(self):
		pass






if __name__=="__main__":
	ut.main()