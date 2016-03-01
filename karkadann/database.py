
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



def get_conn(i_take_responsibility=False):
	if i_take_responsibility:
		#danger danger high voltage
		return mysql.connect(host=host,user=username,passwd=password,db=dbname)
		#This is meant to allow threaded access without having to do with closing(conn)
		#except when you really need it
	if get_conn.conn:
		return get_conn.conn
	else:
		get_conn.conn = mysql.connect(host=host,user=username,passwd=password,db=dbname)
		return get_conn.conn
get_conn.conn = None


