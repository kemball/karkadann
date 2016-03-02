
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
		if not os.path.exists(os.path.join(gb_location,str(salt))):
			SeqIO.write(record,os.path.join(gb_location,str(salt)),'genbank')
			return str(salt)
		else:
			return save_record(record,str(salt)+sample(ascii_lowercase,1))
	else:
		return save_record(record,sample(ascii_lowercaseq,10))

def read_record(assem_id):
	with get_cursor() as cur:
		cur.execute("select gb_record from assemblies where id = %s;",(assem_id,))
		for row in cur:
			return row
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
		return curse.fetchone()



def make_assembly(record,genome_id,reads=None,assembled=None):
	with get_cursor() as curse:
		query = "select id from genomes where id=%s;"
		curse.execute(query,(genome_id,))
		if not curse.fetchall():
			raise ValueError("Tried to make an assembly %s but no matching genome was given")
		gb_location = db.save_record(record)
		query = "insert into assemblies (gb_record,genome_id) values (%s,%s);"
		curse.execute(query,(gb_location,genome_id))
		newassemid = curse.lastrowid
		if reads:
			query = "update assemblies set reads = %s where id = %s;"
			curse.execute(query,(reads,newassemid))
		if assembled:
			query = "update assemblies set assembled = %s where id = %s;"
			curse.execute(query,(assembled,newassemid))





if __name__=="__main__":
	pass
	#tests? who knows


