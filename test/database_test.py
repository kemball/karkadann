from context import karkadann
from karkadann.database import *
import unittest as ut 

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

def record_test():
	#I have no idea how to test if this works, I can't roundtrip with biopython
	pass



def make_assembly_test():
	#wettest code 2016
	try:
		records = SeqIO.parse('/home/kemball/actinobacteria_class/genbank/120971.gb','genbank')
		records= list(records)
		gid = make_genome('test')
	finally:
		with get_cursor() as curse:
			curse.execute("delete from genomes where id=%s;",(gid,))	
	try:
		gid = make_genome('test')
		newid = make_assembly(records,gid)

	finally:
		with get_cursor() as curse:
			curse.execute("delete from assemblies where id = %s",(newid,))
			curse.execute("delete from genomes where id=%s;",(gid,))	

	try:
		gid = make_genome('test')
		newid = make_assembly(records,gid)
		with get_cursor() as curse:
			curse.execute("select gb_record from assemblies where id = %s",(newid,))
			salt = curse.fetchone()[0]
			assert(os.path.exists(os.path.join(gb_location,salt)))
			curse.execute("select genomes.id,genomes.name from genomes\
											 join assemblies on \
											 assemblies.genome_id=genomes.id\
											  where assemblies.id=%s",(newid,))
			fetched_id,fetched_name = curse.fetchone()
			assert(fetched_id==gid)
			assert(fetched_name=='test')
			roundaboutrecord = list(read_record(newid))
			assert(len(roundaboutrecord)==len(records))
			assert(roundaboutrecord[0].id == records[0].id)
			os.remove(os.path.join(gb_location,salt))

	finally:
		with get_cursor() as curse:
			curse.execute("delete from assemblies where id = %s",(newid,))
			curse.execute("delete from genomes where id=%s;",(gid,))	
oldtest = ut.FunctionTestCase(make_assembly_test)


if __name__=="__main__":
	ut.main()