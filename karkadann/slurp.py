


from database import get_cursor as gc 

from Bio import SeqIO
import os
from database import data_location as dataloc


def import_genome(genomename):
	with gc() as curse:
		query = r"insert into genomes (name) values (%s);"
		print query
		print query % genomename
		curse.execute(query,(genomename,))
		curse.connection.commit()

if __name__=="__main__":
	import_genome('test')