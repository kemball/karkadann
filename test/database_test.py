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
		newgenome = Genome(genome_name="thisisatestgenome")
		try:
			newgenome.save()
			self.assertTrue(newgenome.is_real())
			with get_cursor() as cur:
				cur.execute("select id from genomes where name = %s;",('thisisatestgenome',))
				self.assertEqual(cur.fetchone()[0],newgenome.is_real())
		finally:
			newgenome.delete()
		self.assertFalse(newgenome.is_real())

	def test_fetch_genome(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			othergenome=Genome.fetch("testgenome")
			self.assertEqual(othergenome.is_real(),newgenome.is_real())
			self.assertEqual(othergenome._name,newgenome._name)
			self.assertEqual(othergenome._name,"testgenome")
		finally:
			newgenome.delete()

	def test_delete(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			othergenome=Genome.fetch("testgenome")
			othergenome.delete()
			self.assertFalse(newgenome.is_real())
		finally:
			newgenome.delete()

	def test_make_by_id(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			gid = newgenome.is_real()
			othergenome = Genome(db_id=gid)
			self.assertEqual(newgenome.added(),othergenome.added())
			self.assertEqual(gid,newgenome.is_real())
		finally:
			newgenome.delete()

	def test_binomial(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			newgenome.binomial("Heffalump jabberwocky subpsp vorpatius")
			self.assertEqual(newgenome.binomial(),"Heffalump jabberwocky subpsp vorpatius")
		finally:
			newgenome.delete()
		with get_cursor() as cur:
			cur.execute("select binomial from genus_species where genome_id=%s;",newgenome._id)
			self.assertEqual(cur.fetchone(),None)





class AssemblyTest(ut.TestCase):
	#this screams for setup and teardown methods
	#and yet here we are
	records = SeqIO.parse(os.path.join(data_location,'test/testassem.gb'),'genbank')
	records = list(records)

	def test_make_assembly(self):
		new_genome = Genome(genome_name='test')
		new_genome.save()
		try:
			newassem = Assembly(self.records,new_genome)
			newassem.save()
			try:
				self.assertTrue(newassem.record())
				
			finally:
				newassem.delete()
			self.assertFalse(newassem.is_real())
		finally:
			new_genome.delete()

	def test_gb_record(self):
		new_genome= Genome(genome_name="foobarbazqux")
		new_genome.save()
		try:
			newassem = Assembly(self.records,new_genome)
			newassem.save()
			try:
				with get_cursor() as cur:
					try:
						cur.execute("select gb_record from assemblies where id = %s;",(newassem.is_real(),))
						name = cur.fetchone()[0]
						self.assertTrue(os.path.exists(os.path.join(gb_location,name)))
					except:
						print cur._last_executed
						print cur.scroll	
						raise
				with get_cursor() as cur:
					cur.execute("select genomes.id,genomes.name from genomes\
																join assemblies \
																	on assemblies.genome_id=genomes.id\
																where assemblies.id=%s;",(newassem.is_real(),))
					fid,fname = cur.fetchone()
					self.assertEqual(new_genome.is_real(),fid)
					self.assertEqual("foobarbazqux",fname)
				roundaboutrecord = list(newassem.record())
				self.assertEqual(len(roundaboutrecord),len(self.records))
				self.assertEqual(roundaboutrecord[0].id , self.records[0].id)
			finally:
				newassem.delete()
		finally:
			new_genome.delete()

	def test_link_genome(self):
		new_genome = Genome(genome_name="foobarbazqux")
		new_genome.save()
		try:
			new_assem = Assembly(self.records,new_genome)
			new_assem.save()
			try:
				two_assem = Assembly.fetch(new_genome.is_real())
				self.assertEqual(two_assem.is_real(),new_assem.is_real())
				two_rec = list(two_assem.record())
				new_rec = list(new_assem.record())
				self.assertEqual(len(two_rec),len(new_rec))
				self.assertEqual(two_rec[-1].id,new_rec[-1].id)
			finally:
				new_assem.delete()
		finally:
			new_genome.delete()


class ContigTest(ut.TestCase):
	records = SeqIO.parse(os.path.join(data_location,'test/testassem.gb'),'genbank')
	records = list(records)
	def setUp(self):
		self.gid = make_genome('test')
		self.aid = make_assembly(self.records,self.gid)
	def tearDown(self):
		delete_record(aid)
		with get_cursor() as cur:
			cur.execute("delete from assemblies where id = %s;",(self.aid,))
			cur.execute("delete from genomes where id = %s;",(self.gid,))

	def test_contig(self):
		contone = save_contig(self.aid,"ATGCA")
		conttwo = save_contig(self.aid,str(self.records[0].seq))
		contthree = save_from_record(self.aid,self.records[0])
		contig_ids = get_contigs(self.aid)
		self.assertIn(contone,contig_ids)
		self.assertIn(conttwo,contig_ids)
		self.assertIn(contthree,contig_ids)
		self.assertEqual(read_contig_seq(contone),"ATGCA")
		self.assertEqual(read_contig_seq(conttwo),read_contig_seq(contthree))
		with get_cursor() as cur:
			cur.executemany("delete from genes where contig = %s;",[(contone,),(conttwo,),(contthree)])
			cur.executemany("delete from contigs where id = %s;",[(contone,),(conttwo,),(contthree)])
		with get_cursor() as cur:
			cur.execute("select id from contigs where assembly_id=%s;",(self.aid,))
			self.assertEqual(cur.fetchone(),None)



class GeneTest(ut.TestCase):
	records = SeqIO.parse(os.path.join(data_location,'test/testassem.gb'),'genbank')
	records = list(records)
	def setUp(self):
		self.gid = make_genome('test')
		self.aid = make_assembly(self.records,self.gid)
		self.contid = save_contig(self.aid,str(self.records[0].seq),self.records[0].id)

	def tearDown(self):
		delete_record(self.aid)
		with get_cursor() as cur:
			cur.execute("delete from contigs where id=%s;",(self.contid,))
			cur.execute("delete from assemblies where id = %s;",(self.aid,))
			cur.execute("delete from genomes where id = %s;",(self.gid,))

	def test_save_genes(self):
		save_genes(self.contid,self.records[0].features)
		self.assertEqual(str(self.records[0].seq),read_contig_seq(self.contid))
		featlist = [feat for feat in self.records[0].features if feat.type=="CDS"]


		with get_cursor() as cur:
			cur.execute("select count(*) from genes where contig = %s;",(self.contid,))
			self.assertEqual(cur.fetchone()[0],len(featlist))


		with get_cursor() as cur:
			cur.execute("delete from genes where contig = %s;",(self.contid,))

	def test_save_gene(self):
		featlist = [feat for feat in self.records[0].features if feat.type=="CDS"]
		feature = featlist[0]
		translation = feature.qualifiers['translation'][0]
		start = feature.location.start
		end = feature.location.end
		strand = str(feature.location.strand)
		gene_id = save_gene(self.contid,translation,start,end,strand)
		with get_cursor() as cur:
			cur.execute("select id from genes where id = %s;",(gene_id,))
			self.assertIsNot(cur.fetchone(),None)
		with get_cursor() as cur:
			cur.execute("delete from genes where id = %s;",(gene_id,))





class HmmTest(ut.TestCase):
	def test_import(self):
		pass






if __name__=="__main__":
	ut.main()