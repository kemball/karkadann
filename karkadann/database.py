
import MySQLdb as mysql 
import ConfigParser,os
from Bio import SeqIO
from pkg_resources import resource_stream
config_stream = resource_stream(__name__,'karkadann.cfg')
config = ConfigParser.SafeConfigParser()
config.readfp(config_stream)


host     = config.get('database','host')
username = config.get('database','username')
password = config.get('database','password')
dbname   = config.get('database','dbname')

data_location = config.get('data','data_location')

if not os.access(data_location,os.R_OK):
	#panic?
	os.mkdir(data_location)

gb_location = os.path.join(data_location,'genbank')

if not os.access(gb_location,os.R_OK):
	#genbank directory not found
	os.mkdir(gb_location)

hmm_location = os.path.join(data_location,'hmm')

if not os.access(hmm_location,os.R_OK):
	os.mkdir(hmm_location)


from contextlib import contextmanager

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
	# AND I can put in pools later.
	try:
		conn = mysql.connect(host=host,user=username,passwd=password,db=dbname)
		curse = conn.cursor()
		yield curse
		conn.commit()
		#if something bad happened don't commit
		# I checked, it works.
	finally:
		curse.close()
		conn.close()

class dbThing:
	def __init__(self,db_id):
		_id = db_id
	def __init__(self):
		_id = None

	def save(self):
		raise NotImplementedError
	def is_real(self):
		raise NotImplementedError
	def delete(self):
		raise NotImplementedError




class Genome(dbThing):

	@classmethod
	def fetch(cls,genome_name):
		with get_cursor() as cur:
			cur.execute("select id from genomes where name = %s;",(genome_name,))
			return Genome(db_id=cur.fetchone()[0])


	def __init__(self,db_id=None,genome_name=None):
		if db_id:
			with get_cursor() as cur:
				cur.execute("select name,added from genomes where id = %s;",(db_id,))
				self._id = db_id
				self._name = cur.fetchone()[0]
			if genome_name and self._name is not genome_name:
				raise KeyError("Genome %s already exists with id %s."%(self._name,self._id))
		else:
			self._name = genome_name
			self._id = None
			self._added = None

	def save(self):
		if not self.is_real():
			with get_cursor() as cur:
				cur.execute("insert into genomes (name) values (%s);",(self._name,))
				self._id = cur.lastrowid

	def delete(self):
		if self.is_real():
			with get_cursor() as cur:
				cur.execute("delete from genus_species where genome_id=%s;",(self._id,))
				cur.execute("delete from genomes where id = %s;",(self._id,))
		

	def is_real(self):
		with get_cursor() as cur:
			cur.execute("select id from genomes where id = %s;",(self._id,))
			for result in cur:
				return result[0]
			return None

	def added(self):
		if self.is_real():
			with get_cursor() as cur:
				cur.execute("select added from genomes where id =%s;",(self._id,))
				for result in cur:
					return result

		with get_cursor() as cur:
			cur.execute("select binomial from genus_species where genome_id=%s;",(self._id,))
			return cur.fetchone()[0]

	def binomial(self,binomial=None):
		if binomial:
			with get_cursor() as cur:
				cur.execute("insert into genus_species (genome_id,binomial) values (%s,%s);",(self._id,binomial))
		else:
			with get_cursor() as cur:
				cur.execute("select binomial from genus_species where genome_id=%s;",(self._id,))
				return cur.fetchone()[0]



class Assembly(dbThing):


	
	def __init__(self,record=None,genome=None,fastq=None,assembled=None,accession=None,db_id=None):
		if db_id:
			with get_cursor() as cur:
				cur.execute("select * from assemblies where id = %s;",(db_id,))
				self._id,self._gb,self._fastq,self._assembled,self._genome_id,self._acc = cur.fetchone()
		else:
			self._gb = Assembly.save_record(record)
			self._genome_id = genome.is_real()
			self._fastq = fastq
			self._assembled = assembled
			self._acc = accession

	@classmethod
	def fetch(cls,genome_id):
		with get_cursor() as cur:
			cur.execute("select id from assemblies where genome_id = %s;",(genome_id,))
			aid = cur.fetchone()[0]
			if aid:
				return Assembly(db_id=aid)


	def save(self):
		with get_cursor() as cur:
			cur.execute("insert into assemblies (gb_record,fastq,assembled,genome_id,accession) values\
									(%s,%s,%s,%s,%s);",(self._gb,self._fastq,self._assembled,self._genome_id,self._acc))
			self._id = cur.lastrowid

	def is_real(self):
		with get_cursor() as cur:
			cur.execute("select id from assemblies where id=%s;",(self._id,))
			for result in cur:
				return result[0]

	def delete(self):
		if self.is_real():
			with get_cursor() as cur:
				cur.execute("delete from assemblies where id =%s;",(self._id))
			os.remove(os.path.join(gb_location,self._gb))

	def record(self):
		if self.is_real():
			return SeqIO.parse(os.path.join(gb_location,self._gb),'genbank')


	
	@staticmethod
	def save_record(record,salt=None):
		from random import sample
		from string import ascii_lowercase
		if salt:
			if not os.path.exists(os.path.join(gb_location,salt)):
				SeqIO.write(record,os.path.join(gb_location,salt),'genbank')
				return salt
			else:
				return Assembly.save_record(record,salt+str(sample(ascii_lowercase,1)))
		else:
			return Assembly.save_record(record,"".join(sample(ascii_lowercase,10)))






def make_assembly(record,genome_id,reads=None,assembled=None,accession=None):
	new_assembly = Assembly(record,genome=genome_id,fastq=reads,assembled=assembled,accession=accession)
	return new_assembly.is_real()



def save_contig(assembly_id,sequence,accession=None):
	with get_cursor() as cur:
		if accession:
			cur.execute("insert into contigs \
						(assembly_id,sequence,accession)\
						 values (%s,%s,%s);",
						 (assembly_id,sequence,accession))
		else:
			cur.execute("insert into contigs \
				(assembly_id,sequence) values (%s,%s);",(assembly_id,sequence))
		return cur.lastrowid

def save_from_record(assembly_id,record):
	#if record.id isn't the accession number some stuff is going to get a little bit broken.
	contid = save_contig(assembly_id,str(record.seq),record.id)
	save_genes(contid,record.features)
	return contid

def get_contigs(assemid):
	with get_cursor() as curse:
		curse.execute("select id from contigs where assembly_id=%s;",(assemid,))
		return [contig_id[0] for contig_id in curse]


def read_contig_seq(contig_id):
	with get_cursor() as cur:
		query = "select sequence,accession from contigs where id = %s;"
		cur.execute(query,(contig_id,))
		for seq,acc in cur:
			return seq




def save_genes(contig_id,features):
	protlist = [f for f in features if f.type=="CDS"]
	for f in protlist:

		if "protein_id" in f.qualifiers.keys():
			accession = f.qualifiers['protein_id'][0]
		else:
			accession = None

		if "translation" in f.qualifiers.keys():
			translation = f.qualifiers["translation"][0]
		else:
			translation = None
		start = f.location.start
		end = f.location.end
		strand = str(f.location.strand) #str required otherwise the enum in SQL gets confused.
		save_gene(contig_id,translation,start,end,strand,accession)



def save_gene(contig_id,translation,start,end,strand,accession=None):
	with get_cursor() as cur:
		if accession:
			query = "insert into genes (contig,translation,start,end,strand,accession)\
					values (%s,%s,%s,%s,%s,%s);"
			cur.execute(query,(contig_id,translation,start,end,strand,accession))
		else:
			query = "insert into genes (contig,translation,start,end,strand)\
					values (%s,%s,%s,%s,%s);"
			cur.execute(query,(contig_id,translation,start,end,strand))

		return cur.lastrowid


def import_hmm(hmmfile):
	shortname = os.path.basename(hmmfile)
	with get_cursor() as cur:
		query = "insert into hmm_storage filename value %s;"
		cur.execute(query,(shortname,))
	from shutils import copyfile
	copyfile(hmmfile,os.path.join(hmm_location,shortname))
	return shortname


