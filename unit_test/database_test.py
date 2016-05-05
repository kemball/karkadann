import unittest as ut

from karkadann.database import *


class FileTest(ut.TestCase):
	def test_data_location(self):
		self.assertTrue(os.access(data_location, os.R_OK))

	def test_gb_location(self):
		self.assertTrue(os.access(gb_location, os.R_OK))

	def test_hmm_location(self):
		self.assertTrue(os.access(hmm_location, os.R_OK))


class GenomeTest(ut.TestCase):
	def test_make_genome(self):
		newgenome = Genome(genome_name="thisisatestgenome")
		try:
			newgenome.save()
			self.assertTrue(newgenome.is_real())
			with get_cursor() as cur:
				cur.execute("select id from genomes where name = %s;", ('thisisatestgenome',))
				self.assertEqual(cur.fetchone()[0], newgenome.is_real())
		finally:
			newgenome.delete()
		self.assertFalse(newgenome.is_real())

	def test_make_insane_genome(self):

		with self.assertRaises(Exception):
			Genome.get(42)

	def test_fetch_genome(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			othergenome = Genome.fetch("testgenome")
			self.assertEqual(othergenome.is_real(), newgenome.is_real())
			self.assertEqual(othergenome._name, newgenome._name)
			self.assertEqual(othergenome._name, "testgenome")
		finally:
			newgenome.delete()

	def test_delete(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			othergenome = Genome.fetch("testgenome")
			othergenome.delete()
			self.assertFalse(newgenome.is_real())
		finally:
			newgenome.delete()

	def test_make_by_id(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			gid = newgenome.is_real()
			othergenome = Genome.get(gid)
			self.assertEqual(newgenome.added(), othergenome.added())
			self.assertEqual(gid, newgenome.is_real())
		finally:
			newgenome.delete()

	def test_binomial(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			newgenome.binomial("Heffalump jabberwocky subpsp vorpatius")
			self.assertEqual(newgenome.binomial(), "Heffalump jabberwocky subpsp vorpatius")
		finally:
			newgenome.delete()
		with get_cursor() as cur:
			cur.execute("select binomial from genus_species where genome_id=%s;", newgenome._id)
			self.assertEqual(cur.fetchone(), None)

	def test_assemblies(self):
		newgenome = Genome(genome_name="testgenome")
		newgenome.save()
		try:
			self.assertRaises(StopIteration, next, newgenome.assemblies())
		finally:
			newgenome.delete()


class AssemblyTest(ut.TestCase):
	# this screams for setup and teardown methods
	# and yet here we are
	records = SeqIO.parse(os.path.join(data_location, 'test/testassem.gb'), 'genbank')
	records = list(records)

	def test_make_assembly(self):
		new_genome = Genome(genome_name='test_genome')
		new_genome.save()
		try:
			newassem = Assembly(self.records, new_genome)
			newassem.save()
			try:
				self.assertTrue(newassem.record())

			finally:
				newassem.delete()
			self.assertFalse(newassem.is_real())
		finally:
			new_genome.delete()

	def test_parent_watches(self):
		new_genome = Genome(genome_name='test_genome')
		new_genome.save()
		try:
			newassem = Assembly(self.records, new_genome)
			newassem.save()
			try:
				assembly = next(new_genome.assemblies())
				self.assertEqual(assembly.is_real(), newassem.is_real())

			finally:
				newassem.delete()
		finally:
			new_genome.delete()

	def test_delete_by_genome(self):
		new_genome = Genome(genome_name="squiffle")
		new_genome.save()
		new_assem = Assembly(self.records, new_genome)
		new_assem.save()
		new_genome.delete()
		self.assertFalse(new_assem.is_real())

	def test_gb_record(self):
		new_genome = Genome(genome_name="foobarbazqux")
		new_genome.save()
		try:
			newassem = Assembly(self.records, new_genome)
			newassem.save()
			try:
				with get_cursor() as cur:
					try:
						cur.execute("select gb_record from assemblies where id = %s;", (newassem.is_real(),))
						name = cur.fetchone()[0]
						self.assertTrue(os.path.exists(os.path.join(gb_location, name)))
					except:
						print cur._last_executed
						print cur.scroll
						raise
				with get_cursor() as cur:
					cur.execute("select genomes.id,genomes.name from genomes\
																join assemblies \
																	on assemblies.genome_id=genomes.id\
																where assemblies.id=%s;", (newassem.is_real(),))
					fid, fname = cur.fetchone()
					self.assertEqual(new_genome.is_real(), fid)
					self.assertEqual("foobarbazqux", fname)
				roundaboutrecord = list(newassem.record())
				self.assertEqual(len(roundaboutrecord), len(self.records))
				self.assertEqual(roundaboutrecord[0].id, self.records[0].id)
			finally:
				newassem.delete()
		finally:
			new_genome.delete()

	def test_link_genome(self):
		new_genome = Genome(genome_name="foobarbazqux")
		new_genome.save()
		try:
			new_assem = Assembly(self.records, new_genome)
			new_assem.save()
			try:
				two_assem = Assembly.fetch(new_genome.is_real())[0]
				self.assertEqual(two_assem.is_real(), new_assem.is_real())
				two_rec = list(two_assem.record())
				new_rec = list(new_assem.record())
				self.assertEqual(len(two_rec), len(new_rec))
				self.assertEqual(two_rec[-1].id, new_rec[-1].id)
			finally:
				new_assem.delete()
		finally:
			new_genome.delete()


class ContigTest(ut.TestCase):
	records = SeqIO.parse(os.path.join(data_location, 'test/testassem.gb'), 'genbank')
	records = list(records)

	def setUp(self):
		self.test_genome = Genome(genome_name="testgenome")
		self.test_genome.save()
		self.test_assembly = Assembly(record=self.records, genome=self.test_genome)
		self.test_assembly.save()

	def tearDown(self):
		self.test_genome.delete()
		self.assertFalse(self.test_assembly.is_real())

	def test_make_contig(self):
		test_contig = Contig(seq="ATCGATCG", assembly=self.test_assembly)
		test_contig.save()
		try:
			self.assertEqual(test_contig.seq(), "ATCGATCG")
			self.assertTrue(test_contig.is_real())
		finally:
			test_contig.delete()
		self.assertFalse(test_contig.is_real())

	def test_delete_assembly(self):
		other_assem = Assembly(record=self.records, genome=self.test_genome)
		other_assem.save()
		test_contig = Contig(seq="ATCGATCG", assembly=other_assem)
		test_contig.save()
		self.assertTrue(test_contig.is_real())
		other_assem.delete()
		self.assertFalse(test_contig.is_real())

	def test_long_contig(self):
		# about 8Mb
		sequence = "ATCGATCG" * 1000000
		test_contig = Contig(seq=sequence, assembly=self.test_assembly)
		test_contig.save()
		try:
			self.assertTrue(test_contig.seq(), sequence)
		finally:
			test_contig.delete()

	def test_set_accession(self):
		test_contig = Contig(seq="ATCGATCG", assembly=self.test_assembly, accession="WP_1337.1")
		test_contig.save()
		try:
			self.assertEqual(test_contig.acc(), "WP_1337.1")
		finally:
			test_contig.delete()

	def test_get_genes(self):
		test_contig = Contig(seq="ATCGATCG", assembly=self.test_assembly, accession="WP_1337.1")
		test_contig.save()
		try:
			# I know the laziest but genes need their own tests too.
			test_gene = Gene(translation="MICHAELBIOLOGIST",
			                 contig=test_contig,
			                 start=0,
			                 end=2,
			                 strand="-1",
			                 accession="WP_1337.1")
			test_gene.save()
			genes = list(test_contig.genes())
			self.assertEqual(genes[0].is_real(), test_gene.is_real())
			self.assertEqual(len(genes), 1)
		finally:
			test_contig.delete()
		self.assertFalse(test_gene.is_real())


class GeneTest(ut.TestCase):
	# TODO save_many unit_test
	records = SeqIO.parse(os.path.join(data_location, 'test/testassem.gb'), 'genbank')
	records = list(records)

	def setUp(self):
		self.test_genome = Genome(genome_name="testgenome")
		self.test_genome.save()
		self.test_assembly = Assembly(record=self.records, genome=self.test_genome)
		self.test_assembly.save()
		self.test_contig = Contig(seq="ATCGATCG" * 1000, assembly=self.test_assembly)
		self.test_contig.save()

	def tearDown(self):
		self.test_genome.delete()
		self.assertFalse(self.test_assembly.is_real())
		self.assertFalse(self.test_contig.is_real())

	def test_make_gene(self):
		test_gene = Gene(translation="MICHAELBIOLOGIST",
		                 contig=self.test_contig,
		                 start=0,
		                 end=2,
		                 strand=-1,
		                 accession="WP_1337.1")
		test_gene.save()
		try:
			self.assertTrue(test_gene.is_real())
			self.assertEqual(test_gene.location.start, 0)
			self.assertEqual(test_gene.location.strand, -1)
			self.assertEqual(test_gene.translation, "MICHAELBIOLOGIST")
		finally:
			test_gene.delete()

	def test_load_gene(self):
		test_gene = Gene(translation="MICHAELBIOLOGIST",
		                 contig=self.test_contig,
		                 start=0,
		                 end=2,
		                 strand=-1,
		                 accession="WP_1337.1")
		test_gene.save()
		other_gene = Gene.get(test_gene.is_real())
		try:
			self.assertEqual(test_gene.is_real(), other_gene.is_real())
			self.assertEqual(test_gene.__dict__, other_gene.__dict__)
		finally:
			test_gene.delete()

	def test_record(self):
		test_gene = Gene(translation="MICHAELBIOLOGIST",
		                 contig=self.test_contig,
		                 start=0,
		                 end=2,
		                 strand=-1,
		                 accession="WP_1337.1")
		test_gene.save()
		SeqIO.write(test_gene.record, 'foo.fasta', 'fasta')
		round_rec = SeqIO.read('foo.fasta', 'fasta')
		self.assertEqual(round_rec.id, test_gene.record.id)
		self.assertEqual(round_rec.seq, test_gene.record.seq)
		os.remove('foo.fasta')


class HitTest(ut.TestCase):
	# TODO save_many unit_test

	from karkadann.assimilate import assimilate_from_ncbi

	testgbfile = os.path.join(data_location, 'test/testassem.gb')
	ng = assimilate_from_ncbi(testgbfile)
	assem = next(ng.assemblies())
	contig = next(assem.contigs())

	@classmethod
	def tearDownClass(cls):
		cls.ng.delete()

	def test_make_hit(self):
		gene = next(self.contig.genes())
		newhit = Hit(gene=gene, score=5, hmm="testhmm", seq="ACDEFGHIKLMNPQRSTVWY")
		try:
			newhit.save()
			self.assertTrue(newhit.is_real())
			self.assertEqual(5, newhit.score)
			otherhit = Hit.get(newhit.is_real())
			self.assertEqual(otherhit.score, newhit.score)
			otherhit.save()
			self.assertEqual(otherhit.is_real(), newhit.is_real())
			self.assertEqual(str(newhit._seq), "ACDEFGHIKLMNPQRSTVWY")
		finally:
			newhit.delete()

	def test_fetch_hit(self):
		gene = next(self.contig.genes())
		newhit = Hit(gene=gene, score=5, hmm="testhmm")
		try:
			newhit.save()
			otherhit = Hit.fetch(gene)
			self.assertEqual(otherhit[0].is_real(), newhit.is_real())
		finally:
			newhit.delete()


class ClusterTest(ut.TestCase):
	# TODO save_many test

	from karkadann.assimilate import assimilate_from_ncbi
	testgbfile = os.path.join(data_location, 'test/testassem.gb')
	ng = assimilate_from_ncbi(testgbfile)
	assem = next(ng.assemblies())
	contig = next(assem.contigs())

	@classmethod
	def tearDownClass(cls):
		cls.ng.delete()

	def test_make_cluster(self):
		newcluster = Cluster(gene_list=list(self.contig.genes()), classification="testclass")
		try:
			newcluster.save()
			self.assertTrue(newcluster.is_real())
			othercluster = Cluster.get(db_id=newcluster.is_real())
			self.assertTrue(othercluster.is_real())
			self.assertEqual(newcluster._kind, othercluster._kind)
			self.assertEqual(newcluster._kind, 'testclass')
		finally:
			newcluster.delete()
		self.assertFalse(newcluster.is_real())

	def test_get_fna(self):
		newcluster = Cluster(gene_list=list(self.contig.genes()), classification="testclass")
		newcluster.save()
		try:
			records = newcluster.fna()
			from tempfile import NamedTemporaryFile as ntf
			with ntf(mode='w', delete=True) as tempfile:
				SeqIO.write(records, tempfile.name, 'fasta')
				roundtrip = SeqIO.parse(tempfile.name, 'fasta')
				roundtrip = list(roundtrip)[0]
				self.assertEqual(len(roundtrip), len(records))
				self.assertEqual(roundtrip.id, records.id)
		finally:
			newcluster.delete()

	def test_get_faa(self):
		glist = list(self.contig.genes())
		newcluster = Cluster(gene_list=glist, classification="testclass")
		newcluster.save()
		try:
			records = newcluster.faa()
			othercluster = Cluster.get(db_id=newcluster.is_real())
			otherrecords = othercluster.faa()
			self.assertEqual(othercluster._kind, newcluster._kind)
			self.assertEqual(othercluster._kind, 'testclass')
			self.assertEqual(newcluster.is_real(), othercluster.is_real())
			self.assertEqual(len(records), len(otherrecords))
			self.assertEqual(str(records[0]), str(otherrecords[0]))
		finally:
			newcluster.delete()
		self.assertFalse(newcluster.is_real())
		self.assertFalse(othercluster.is_real())


if __name__ == "__main__":
	ut.main()
