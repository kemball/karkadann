


from database import get_cursor as gc 
from database import data_location as dataloc
import database as db

from Bio import SeqIO
import os

def slurp_genbank(gbfilename):
	pass

def slurp_fasta(fastafilename):
	pass





if __name__=="__main__":
	genome_id = db.make_genome('test')
	print db.get_genome("squiffle")
	print db.get_genome("test")
