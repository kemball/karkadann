
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
gb_location = os.path.join(data_location,'genbank')

if not os.access(data_location,os.R_OK):
	#panic?
	os.mkdir(data_location)

if not os.access(gb_location,os.R_OK):
	#genbank directory not found
	os.mkdir(gb_location)


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
		# I haven't actually tested whether that works. :(
	finally:
		curse.close()
		conn.close()

from random import sample
from string import ascii_lowercase
def save_record(record,salt=None):
	if salt:
		if not os.path.exists(os.path.join(gb_location,salt)):
			SeqIO.write(record,os.path.join(gb_location,salt),'genbank')
			return salt
		else:
			return save_record(record,salt+str(sample(ascii_lowercase,1)))
	else:
		return save_record(record,"".join(sample(ascii_lowercase,10)))

def read_record(assem_id):
	with get_cursor() as cur:
		cur.execute("select gb_record from assemblies where id = %s;",(assem_id,))
		for salt in cur:
			return SeqIO.parse(os.path.join(gb_location,salt[0]),'genbank')
		raise IOError("assembly %s has no gb_record"%assem_id)


def make_genome(genomename):
	with get_cursor() as curse:
		query = r"insert into genomes (name) values (%s);"
		curse.execute(query,(genomename,))
		return curse.lastrowid

def get_genome(genomename):
	with get_cursor() as curse:
		query = "select id from genomes where name = %s;"
		curse.execute(query,(genomename,))
		result = curse.fetchall()
		for q in result:
			return q[0]
		return None

def genome_test():
	with get_cursor() as cursorone:
		with get_cursor() as cursortwo:
			assert( id(cursorone) is not id(cursortwo))

	testid = make_genome("thisisnotarealname")
	assert(testid)
	with get_cursor() as cur:
		cur.execute("delete from genomes where id = %s;",(testid,))
	testid = make_genome("thisisnotarealname")
	retrievedid = get_genome("thisisnotarealname")
	assert(testid==retrievedid)
	with get_cursor() as cur:
		cur.execute("delete from genomes where id = %s",(testid,))
		#so that shouldn't commit until it falls otu of scope, right?
		assert( get_genome("thisisnotarealname"))
	with get_cursor() as cursorone:
		cursorone.execute("insert into genomes (name,id) values (%s,%s);",("thisisnotarealname",42))
		with get_cursor() as cursortwo:
			assert(not get_genome("thisisnotarealname"))






def make_assembly(record,genome_id,reads=None,assembled=None,accession=None):
	with get_cursor() as curse:
		query = "select id from genomes where id=%s;"
		curse.execute(query,(genome_id,))
		if not curse.fetchall():
			raise ValueError("Tried to make an assembly %s but no matching genome was given")
		gb_location = save_record(record)
		query = "insert into assemblies (gb_record,genome_id) values (%s,%s);"
		curse.execute(query,(gb_location,genome_id))
		newassemid = curse.lastrowid
		if reads:
			query = "update assemblies set reads = %s where id = %s;"
			curse.execute(query,(reads,newassemid))
		if assembled:
			query = "update assemblies set assembled = %s where id = %s;"
			curse.execute(query,(assembled,newassemid))
		if accession:
			query = "update assemblies set accession = %s where id = %s"
			curse.execute(query,(accession,newassemid))
		return curse.lastrowid

def make_assembly_test():
	try:
		records = SeqIO.parse('/home/kemball/actinobacteria_class/genbank/120971.gb','genbank')
		records= list(records)
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




if __name__=="__main__":
	genome_test()

	make_assembly_test()
	#tests? who knows


