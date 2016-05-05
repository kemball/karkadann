import ConfigParser
import MySQLdb as mysql
import os

from Bio.Alphabet import IUPAC

from Bio import SeqFeature
from Bio import SeqIO, SeqRecord
from Bio.Seq import Seq
from contextlib import contextmanager

from Bio.Seq import Seq
from pkg_resources import resource_stream

config_stream = resource_stream(__name__, 'karkadann.cfg')
config = ConfigParser.SafeConfigParser()
config.readfp(config_stream)

host = config.get('database', 'host')
username = config.get('database', 'username')
password = config.get('database', 'password')
dbname = config.get('database', 'dbname')

data_location = config.get('data', 'data_location')

if not os.access(data_location, os.R_OK):
	os.mkdir(data_location)

gb_location = os.path.join(data_location, 'genbank')

if not os.access(gb_location, os.R_OK):
	# genbank directory not found
	os.mkdir(gb_location)

hmm_location = os.path.join(data_location, 'hmm')

if not os.access(hmm_location, os.R_OK):
	os.mkdir(hmm_location)

import multiprocessing
# TODO: a config option for this.
try:
	num_cores = multiprocessing.cpu_count()
except NotImplementedError:
	num_cores = 5

@contextmanager
def get_cursor():
	"""returns a cursor that can have .execute and .fetchall and so on called on it in a thread-safe way"""
	# That thread-safety bears some explaining. because of contextmanager magic
	# get_cursor() returns a generator with a single item in it: the connection
	# The generator wrapper lets us call .close() appropriately when the generator returns
	# this is for use in `with get_cursor() as curse:`

	conn = mysql.connect(host=host, user=username, passwd=password, db=dbname)
	curse = conn.cursor()
	try:
		yield curse
		conn.commit()

	# if something bad happened don't commit
	# I checked, it works.
	finally:
		curse.close()
		conn.close()
	# this makes a new connect *every* time.
	# it does goddamn work tho


def cursor_required(f):
	# this also bears some explaining. Quite a few of the .save methods
	# and suchlike need a cursor to the database, but it's pretty inefficient
	# to open and close it for every interaction. The methods need to be called
	# without keyword arguments so other modules can use them without touching
	# get_cursor(), but I didn't want to write these three lines in every
	# single method. So I wrote a decorator to do it for me.
	# This way database actions can share a cursor to their dependents
	# but it's not compulsory.
	def cursored_fxn(*args, **kwargs):
		# enforce keyworded cur argument.
		if 'cur' in kwargs.keys():
			return f(*args, **kwargs)
		else:
			with get_cursor() as cur:
				return f(*args, cur=cur, **kwargs)

	return cursored_fxn


class dbThing(object):
	def save(self, cur=None):
		raise NotImplementedError

	def is_real(self, cur=None):
		raise NotImplementedError

	def delete(self, cur=None):
		raise NotImplementedError

	@classmethod
	def get(cls, db_id):
		return cls.get_many([db_id]).next()

	@classmethod
	def get_many(cls, db_ids):
		raise NotImplementedError


class Genome(dbThing):
	@classmethod
	def fetch(cls, genome_name):
		with get_cursor() as cur:
			cur.execute("select id from genomes where name = %s;", (genome_name,))
			return Genome.get(cur.fetchone()[0])

	@classmethod
	def get_many(cls, db_ids):
		with get_cursor() as cur:
			for db_id in db_ids:
				cur.execute("select id,name,added from genomes where id = %s;", (db_id,))
				gid, name, added = cur.fetchone()
				ng = Genome()
				ng._id,ng._name, ng._added = gid,name, added
				yield ng

	def __init__(self, genome_name=None):
		self._name = genome_name
		self._id = None
		self._added = None

	@cursor_required
	def save(self, cur=None):
		if not self._id:
			cur.execute("insert into genomes (name) values (%s);", (self._name,))
			self._id = cur.lastrowid

	@cursor_required
	def delete(self, cur=None):
		if self.is_real(cur=cur):
			for assem in self.assemblies():
				assem.delete(cur=cur)
			cur.execute("delete from genomes where id = %s;", (self._id,))

	@cursor_required
	def is_real(self, cur=None):
		if self._id:
			cur.execute("select id from genomes where id = %s;", (self._id,))
			for result in cur.fetchall():
				return result[0]
		return None

	@cursor_required
	def added(self, cur=None):
		if self.is_real(cur=cur):
			cur.execute("select added from genomes where id =%s;", (self._id,))
			for result in cur.fetchall():
				return result

	@cursor_required
	def binomial(self, binomial=None, cur=None, ):
		if binomial:
			cur.execute("insert into genus_species (genome_id,binomial) values (%s,%s);", (self._id, binomial))
		else:
			cur.execute("select binomial from genus_species where genome_id=%s;", (self._id,))
			return cur.fetchone()[0]


	def assemblies(self):
		with get_cursor() as cur:
			if self.is_real(cur=cur):
				cur.execute("select id from assemblies where genome_id = %s;", (self._id,))
				ids = list([x for (x,) in cur.fetchall()])
				return Assembly.get_many(ids)


class Assembly(dbThing):
	def __init__(self, record=None, genome=None, fastq=None, assembled=None, accession=None):
		self._id = None
		self._gb = record and Assembly.save_record(record) or None
		self._genome_id = genome and genome.is_real() or None
		self._fastq = fastq
		self._assembled = assembled
		self._acc = accession or ""

	@classmethod
	def get_many(cls, db_ids):
		with get_cursor() as cur:
			for db_id in db_ids:
				cur.execute("select * from assemblies where id = %s;", (db_id,))
				nself = cls()
				nself._id, nself._gb, nself._fastq, nself._assembled, nself._genome_id, nself._acc = cur.fetchone()
				yield nself

	@classmethod
	@cursor_required
	def fetch(cls, genome_id, cur=None):
		cur.execute("select id from assemblies where genome_id = %s;", (genome_id,))
		aid = [x for (x,) in cur.fetchall()]
		return list(Assembly.get_many(aid))

	@cursor_required
	def save(self, cur=None):
		if not self._id:
			cur.execute("insert into assemblies (gb_record,fastq,assembled,genome_id,accession) values\
									(%s,%s,%s,%s,%s);",
			            (self._gb, self._fastq, self._assembled, self._genome_id, self._acc))
			self._id = cur.lastrowid

	@cursor_required
	def is_real(self, cur=None):
		if self._id:
			cur.execute("select id from assemblies where id=%s;", (self._id,))
			for result in cur.fetchall():
				return result[0]
		else:
			return False

	@cursor_required
	def delete(self, cur=None):

		if self.is_real(cur=cur):
			cur.execute("delete from assemblies where id =%s;", (self._id,))
			os.remove(os.path.join(gb_location, self._gb))

	@cursor_required
	def record(self, cur=None):

		if self.is_real(cur=cur):
			return SeqIO.parse(os.path.join(gb_location, self._gb), 'genbank')

	@cursor_required
	def contigs(self, cur=None):

		if self.is_real(cur=cur):
			cur.execute("select id from contigs where assembly_id=%s;", (self._id,))
			results = cur.fetchall()
			return (Contig.get(db_id=x[0]) for x in results)

	@staticmethod
	def save_record(record, salt=None):
		from random import sample
		from string import ascii_lowercase
		if salt:
			if not os.path.exists(os.path.join(gb_location, salt)):
				SeqIO.write(record, os.path.join(gb_location, salt), 'genbank')
				return salt
			else:
				return Assembly.save_record(record, salt + str(sample(ascii_lowercase, 1)))
		else:
			return Assembly.save_record(record, "".join(sample(ascii_lowercase, 10)))


class Contig(dbThing):
	def __init__(self, seq=None, assembly=None, accession=None):
		self._seq = seq
		self._assembly_id = assembly and assembly.is_real() or None
		self._acc = accession or ""

	@classmethod
	def get_many(cls, db_ids):
		with get_cursor() as cur:
			for db_id in db_ids:
				cur.execute("select id,sequence,assembly_id,accession from contigs where id = %s;", (db_id,))
				ncont = cls()
				ncont._id, ncont._seq, ncont._assembly_id, ncont._acc = cur.fetchone()
				yield ncont

	@cursor_required
	def save(self, cur=None):

		cur.execute("insert into contigs (sequence,assembly_id,accession) \
					  values (%s,%s,%s);", (self._seq, self._assembly_id, self._acc))
		self._id = cur.lastrowid

	@cursor_required
	def is_real(self, cur=None):
		try:
			cur.execute("select id from contigs where id=%s;", (self._id,))
		except AttributeError:
			return False
		for result in cur.fetchall():
			return result[0]

	@cursor_required
	def delete(self, cur=None):

		if self.is_real(cur=cur):
			cur.execute("delete from contigs where id = %s;", (self._id,))

	def seq(self):
		return self._seq

	def acc(self):
		return self._acc

	@cursor_required
	def genes(self, cur=None):
		if self.is_real(cur=cur):
			cur.execute("select id from genes where contig = %s order by start;", (self._id,))
			return Gene.get_many([x for (x,) in cur.fetchall()])


class Gene(dbThing):
	def __init__(self,
	             translation=None,
	             contig=None,
	             start=None,
	             end=None,
	             strand=None,
	             accession=None):
		self._id = None
		self._translation = translation
		self._contig = contig and contig.is_real() or None
		self._start = start
		self._end = end
		self._strand = strand
		self._acc = accession or ""

	@classmethod
	def get_many(cls, db_ids):
		if not len(db_ids):
			return
		with get_cursor() as cur:
			query = "select id,translation,contig,start,end,strand,accession from genes where id in (%s);"
			query %= ",".join(['%s']*len(db_ids))
			try:
				cur.execute(query,tuple(db_ids))
			except:
				print "something went wrong??"
				print query%tuple(db_ids)
				raise
			for res in cur.fetchall():
				ng = Gene()
				ng._id, ng._translation, ng._contig, ng._start, ng._end, ng._strand, ng._acc = res
				ng._start = int(ng._start)
				ng._end = int(ng._end)
				ng._strand = int(ng._strand)
				yield ng

	@cursor_required
	def save(self, cur=None):
		if not self._id:
			cur.execute("insert into genes (translation,start,end,strand,contig,accession) values \
			                            (%s,%s,%s,%s,%s,%s);",
			            (self._translation, self._start, self._end, self._strand, self._contig, self._acc))
			self._id = cur.lastrowid

	@staticmethod
	def _save_many(genes_iterable):
		with get_cursor() as cur:
			tuples = [(x._translation, x._start, x._end, x._strand, x._contig, x._acc) for x in genes_iterable]
			cur.executemany(
				"insert into genes (translation, start, end, strand, contig, accession) values (%s,%s,%s,%s,%s,%s);",
				tuples)

	@cursor_required
	def is_real(self, cur=None):
		if self._id:
			cur.execute("select id from genes where id = %s;", (self._id,))
			for result in cur.fetchall():
				return result[0]
		else:
			return False

	@cursor_required
	def delete(self, cur=None):
		if self.is_real(cur=cur):
			cur.execute("delete from genes where id =%s;", (self._id,))

	@property
	def translation(self):
		return self._translation

	@translation.setter
	def translation(self, new_trans):
		self._translation = str(new_trans)

	@property
	def location(self):
		return SeqFeature.FeatureLocation(int(self._start), int(self._end), strand=int(self._strand))

	@location.setter
	def location(self, new_location):
		self._start = new_location.start
		self._end = new_location.end
		self._strand = new_location.strand

	@property
	def record(self):
		return SeqRecord.SeqRecord(id="%s|%s_%s" % (self._contig, self._contig, self._id),
		                           seq=Seq(self._translation, alphabet=IUPAC.protein))


	def hits(self):
		with get_cursor() as cur:
			cur.execute("select id from hits where gene = %s;",(self.is_real(),))
			return list(Hit.get_many([x for (x,) in cur.fetchall()]))

	def orthogroup(self, batch=None):
		if self.is_real():
			with get_cursor() as cur:
				if not batch:
					batch = most_recent_batch(cur=cur)
				cur.execute("select `group_name` from orthogroups where gene =%s and batch=%s;", (self._id, batch))
				for result in cur:
					return result[0]
		else:
			return ""


class Hit(dbThing):
	def __init__(self, gene=None, score=None, seq=None, hmm=None):
		self._seq = seq and Seq(seq,IUPAC.protein)
		self._score = score
		self._hmm = hmm
		self._gene = gene and gene.is_real() or None
		self._id = None

	@classmethod
	def get_many(cls, db_ids):
		if not db_ids:
			return
		with get_cursor() as cur:
			query = "select id,score,gene,hmm,hitseq from hits where id in (%s);"
			query %= ",".join(['%s']*len(db_ids))
			try:
				cur.execute(query,tuple(db_ids))
			except:
				print query%tuple(db_ids)
				print "PANIC SOMETHING BAD HAPPENED"
				raise
			for (ID,SCORE,GENE,HMM,SEQ) in cur.fetchall():
				nhit = Hit(score=SCORE, seq=SEQ,hmm=HMM)
				# An additional dependence?
				# TODO think a little bit about whether ._gene and the like should be integers or Genes.
				nhit._gene = GENE
				nhit._id= ID
				yield nhit

	@classmethod
	def fetch(cls, gene):
		with get_cursor() as cur:
			cur.execute("select id from hits where gene = %s;", (gene.is_real(),))
			return list(Hit.get_many([x for (x,) in cur.fetchall()]))

	@cursor_required
	def delete(self, cur=None):
		cur.execute("delete from hits where id = %s;", (self._id,))

	@cursor_required
	def save(self, cur=None):
		if not self.is_real():
			cur.execute("insert into hits (score,gene,hmm,hitseq) values (%s,%s,%s,%s);",
			            (self._score, self._gene, self._hmm, self._seq))
			self._id = cur.lastrowid

	@staticmethod
	def _save_many(hits_iterable):
		with get_cursor() as cur:
			tuples = [(x._score, x._gene, x._hmm, x._seq) for x in hits_iterable]
			cur.executemany("insert into hits (score,gene,hmm,hitseq) values (%s,%s,%s,%s);", tuples)

	@cursor_required
	def is_real(self, cur=None):
		if self._id:
			cur.execute("select id from hits where id = %s;", (self._id,))
			for result in cur:
				return result[0]
			return None

	@property
	def score(self):
		return self._score

	@property
	def seq(self):
		return self._seq

	@property
	def hmm(self):
		return self._hmm


class Cluster(dbThing):
	def __init__(self, gene_list=[], classification=None):

		self._gene_list = gene_list
		self._id = None
		self._kind = classification

	@classmethod
	def get_many(cls, db_ids):
		with get_cursor() as cur:
			for db_id in db_ids:
				ncl = Cluster()
				cur.execute("select gene,classification from clusters where id = %s;", (db_id,))
				gids = []
				for g, c in cur.fetchall():
						gids.append(g)
						if c:
							ncl._kind = c
				ncl._gene_list = list(Gene.get_many(gids))
				ncl._id = db_id
				yield ncl

	@property
	def gene_list(self):
		return self._gene_list

	@cursor_required
	def is_real(self, cur=None):
		if self._id:
			cur.execute("select id from cluster_names where id = %s;", (self._id,))
			for r in cur.fetchall():
				return self._id
		return None

	@cursor_required
	def save(self, cur=None):
		if not self.is_real():
			# doesn't have to be unique
			cur.execute("select genomes.id from genomes join assemblies on assemblies.genome_id=genomes.id "
			            "join contigs on contigs.assembly_id=assemblies.id "
			            "join genes on genes.contig=contigs.id where genes.id=%s;", (self._gene_list[0].is_real(),))
			genome_id = cur.fetchone()[0]
			# TODO munge this for uniqueness and QOL
			cur.execute("insert into cluster_names (name,genome) values (%s,%s);", (self._kind, genome_id))
			self._id = cur.lastrowid
			for gene in self._gene_list:
				cur.execute("insert into clusters (id,classification,gene) values (%s,%s,%s);",
				            (self._id, self._kind, gene.is_real()))

	@cursor_required
	def delete(self, cur=None):
		if self.is_real():
			cur.execute("delete from cluster_names where id = %s;", (self._id,))

	@cursor_required
	def faa(self, cur=None):
		if self.is_real(cur=cur):
			cur.execute("select genes.id,genes.translation from "
			            "genes join clusters on "
			            "genes.id=clusters.gene where clusters.id=%s order by genes.start;",
			            (self._id,))
			return [SeqRecord.SeqRecord(id=str(gene_id), seq=Seq(trans, alphabet=IUPAC.protein)) for gene_id, trans in
			        cur.fetchall()]

	@cursor_required
	def fna(self, cur=None):
		if self.is_real(cur=cur):
			beginning = min([gene.location.start for gene in self._gene_list])
			end = max([gene.location.end for gene in self._gene_list])
			cur.execute("select sequence from contigs join genes on genes.contig =contigs.id where genes.id=%s;",
			            (self._gene_list[0].is_real(),))
			contigseq = Seq(cur.fetchone()[0], alphabet=IUPAC.IUPACAmbiguousDNA)[beginning:end]
			contig = SeqRecord.SeqRecord(seq=contigseq, id=str(self._id))
			return contig


@cursor_required
def most_recent_batch(cur=None):
	cur.execute("select id from orthomcl_batches order by done desc limit 1;")
	for result in cur:
		return result[0]


@cursor_required
def start_batch(cur=None):
	# this looks weird but the defaults in the table are what I actually want.
	# go figure.
	cur.execute("insert into orthomcl_batches () values ();")
	return cur.lastrowid


def save_orthogroup(gene, orthogroup, batch=None):
	print "saving orthogroup"
	with get_cursor() as cur:
		if not batch:
			batch = most_recent_batch(cur=cur)
		if not gene.is_real():
			raise ValueError("gene %s needs to be real" % gene._acc)
		cur.execute("insert into orthogroups (batch,`group_name`,gene) values (%s,%s,%s);",
		            (batch, orthogroup, gene.is_real()))
