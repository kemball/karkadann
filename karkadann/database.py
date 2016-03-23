import ConfigParser
import MySQLdb as mysql
import os
import threading
from Bio import SeqFeature
from Bio import SeqIO, SeqRecord
from contextlib import contextmanager
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

threadLocal = threading.local()



@contextmanager
def get_cursor():
	"""returns a cursor that can have .execute and .fetchall and so on called on it in a thread-safe way"""
	# That thread-safety bears some explaining. because of contextmanager magic
	# get_cursor() returns a generator with a single item in it: the connection
	# The generator wrapper lets us call .close() appropriately when the generator returns
	# this is for use in `with get_cursor() as curse:`
	# Generators in general are not thread-safe: if a generator is trying to generate another value
	# after yield returns and is prompted to do that same thing again it gets confused and sad.
	# This looks like a generator, but it's not. It returns many small generators instead of one large
	# connection generator. I could also put connection pools etc in here but that's work so I haven't.
	# the whole point being I don't have to use closing() everywhere, I have thread-safe cursors,
	# AND I can put in pools later. (I already tried once and screwed it up somehow. Oh well.)

	if hasattr(threadLocal, 'conn'):
		conn = threadLocal.conn
	else:
		conn = mysql.connect(host=host, user=username, passwd=password, db=dbname)
		threadLocal.conn = conn
	curse = conn.cursor()
	try:
		yield curse
		conn.commit()
	# if something bad happened don't commit
	# I checked, it works.

	finally:
		curse.close()


class dbThing:
	def save(self, cur=None):
		raise NotImplementedError

	def is_real(self, cur=None):
		raise NotImplementedError

	def delete(self, cur=None):
		raise NotImplementedError


def cursor_required(f):
	def cursored_fxn(*args, **kwargs):
		# enforce keyworded cur argument.
		if 'cur' in kwargs.keys():
			return f(*args, **kwargs)
		else:
			with get_cursor() as cur:
				return f(*args, cur=cur, **kwargs)

	return cursored_fxn


class Genome(dbThing):
	@classmethod
	def fetch(cls, genome_name):
		with get_cursor() as cur:
			cur.execute("select id from genomes where name = %s;", (genome_name,))
			return Genome(db_id=cur.fetchone()[0])

	def __init__(self, db_id=None, genome_name=None):
		if db_id:
			with get_cursor() as cur:
				cur.execute("select name,added from genomes where id = %s;", (db_id,))
				self._id = db_id
				self._name = cur.fetchone()[0]
			if genome_name and self._name is not genome_name:
				raise KeyError("Genome %s already exists with id %s." % (self._name, self._id))
		else:
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
			for assem in self.assemblies(cur=cur):
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

	@cursor_required
	def assemblies(self, cur=None):
		if self.is_real(cur=cur):
			cur.execute("select id from assemblies where genome_id = %s;", (self._id,))
			return (Assembly(db_id=x) for (x,) in cur.fetchall())


class Assembly(dbThing):
	def __init__(self, record=None, genome=None, fastq=None, assembled=None, accession=None, db_id=None):
		if db_id:
			with get_cursor() as cur:
				cur.execute("select * from assemblies where id = %s;", (db_id,))
				self._id, self._gb, self._fastq, self._assembled, self._genome_id, self._acc = cur.fetchone()
		else:
			self._id = None
			self._gb = Assembly.save_record(record)
			self._genome_id = genome.is_real()
			self._fastq = fastq
			self._assembled = assembled
			self._acc = accession or ""

	@classmethod
	@cursor_required
	def fetch(cls, genome_id, cur=None):
		cur.execute("select id from assemblies where genome_id = %s;", (genome_id,))
		aid = cur.fetchall()
		return [Assembly(db_id=x) for (x,) in aid]

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
			cur.execute("delete from assemblies where id =%s;", (self._id))
			os.remove(os.path.join(gb_location, self._gb))

	@cursor_required
	def record(self, cur=None):

		if self.is_real(cur=cur):
			return SeqIO.parse(os.path.join(gb_location, self._gb), 'genbank')

	@cursor_required
	def contigs(self, cur=None):

		if self.is_real(cur=cur):
			cur.execute("select id from contigs where assembly_id=%s;", (self._id,))
			return (Contig(db_id=x[0]) for x in cur.fetchall())

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
	def __init__(self, seq=None, assembly=None, accession=None, db_id=None):
		with get_cursor() as cur:
			if db_id:
				cur.execute("select id,sequence,assembly_id,accession from contigs where id = %s;", (db_id,))
				self._id, self._seq, self._assembly_id, self._acc = cur.fetchone()
			else:
				self._seq = seq
				self._assembly_id = assembly.is_real(cur=cur)
				self._acc = accession or ""

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
			return (Gene(db_id=x) for (x,) in cur.fetchall())


class Gene(dbThing):
	def __init__(self, db_id=None,
	             translation=None,
	             contig=None,
	             start=None,
	             end=None,
	             strand=None,
	             accession=None):
		if db_id:
			with get_cursor() as cur:
				cur.execute("select id,translation,contig,start,end,strand,accession from genes where id=%s;", (db_id,))
				self._id, self._translation, self._contig, self._start, self._end, self._strand, self._acc = cur.fetchone()
				self._start = int(self._start)
				self._end = int(self._end)
				self._strand = int(self._strand)
		else:
			self._id = None
			self._translation = translation
			self._contig = contig.is_real()
			self._start = start
			self._end = end
			self._strand = strand
			self._acc = accession or ""

	@cursor_required
	def save(self, cur=None):
		if not self._id:
			cur.execute("insert into genes (translation,start,end,strand,contig,accession) values \
			                            (%s,%s,%s,%s,%s,%s);",
			            (self._translation, self._start, self._end, self._strand, self._contig, self._acc))
			self._id = cur.lastrowid

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

	def hit_scores(self):
		with get_cursor() as cur:
			cur.execute("select score,hmm from hits where gene = %s;", (self.is_real(),))
			return list(cur.fetchall())


class Hit(dbThing):
	def __init__(self, db_id=None, gene=None, score=None, hmm=None):
		if db_id:
			with get_cursor() as cur:
				cur.execute("select id,score,gene,hmm from hits where id = %s;", (db_id,))
				self._id, self._score, self._gene, self._hmm = cur.fetchone()
		else:
			self._score = score
			self._hmm = hmm
			self._gene = gene.is_real()
			self._id = None

	@classmethod
	def fetch(cls, gene):
		with get_cursor() as cur:
			cur.execute("select id from hits where gene = %s;",(gene.is_real(),))
			return [Hit(db_id=x) for x in cur.fetchall()]

	@cursor_required
	def delete(self, cur=None):
		cur.execute("delete from hits where id = %s;", (self._id,))

	@cursor_required
	def save(self, cur=None):
		if not self.is_real():
			cur.execute("insert into hits (score,gene,hmm) values (%s,%s,%s);", (self._score, self._gene, self._hmm))
			self._id = cur.lastrowid

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


class Cluster(dbThing):

	def __init__(self,db_id=None, gene_list=None, classification=None):
		if db_id:
			self.gene_list = []
			with get_cursor() as cur:
				cur.execute("select gene,classification from clusters where id = %s;", (db_id,))
				for g,c in cur.fetchall():
					self.gene_list.append(Gene(db_id=g))
					if c:
						self._kind = c
				self._id = db_id
		else:
			self.gene_list = gene_list
			self._id = None
			self._kind = classification


	@cursor_required
	def is_real(self, cur=None):
		if self._id:
			cur.execute("select id from clusters where id = %s;", (self._id,))
			for r in cur.fetchall():
				return self._id
		return None

	@cursor_required
	def save(self, cur=None):
		if not self.is_real():
			cur.execute("select count(distinct id) from clusters where classification=%s;", (self._kind,))
			# blah blah race condition blah blah
			self._id = self._kind + 'q' + str(int(cur.fetchone()[0])+1)
			cur.execute("select count(*) from clusters where id=%s;", (self._id,))
			assert(cur.fetchone()[0] == 0)
			for gene in self.gene_list:
				cur.execute("insert into clusters (id,classification,gene) values (%s,%s,%s);",
				                (self._id, self._kind, gene.is_real()))

	@cursor_required
	def delete(self, cur=None):
		if self.is_real():
			cur.execute("delete from clusters where id = %s;", (self._id,))

