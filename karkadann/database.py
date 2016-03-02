
import MySQLdb as mysql 
import ConfigParser,os
from pkg_resources import resource_stream
config_stream = resource_stream(__name__,'karkadann.cfg')
config = ConfigParser.SafeConfigParser()
config.readfp(config_stream)


host     = config.get('database','host')
username = config.get('database','username')
password = config.get('database','password')
dbname   = config.get('database','dbname')

data_location = config.get('data','data_location')
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
	finally:
		curse.close()
		conn.close()

if __name__=="__main__":
	pass
	#tests? who knows


