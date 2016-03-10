


from database import get_cursor as gc 
from database import data_location as dataloc
import database as db

from Bio import SeqIO
import os

def annotate_save(gbfilename,genome_name):
	from prodigal import annotate
	records = list(SeqIO.parse(gbfilename,'genbank'))
	genid = db.make_genome(genome_name)
	organisms = set([])
	for r in records:
		try:
			organisms |= set([r.annotations['organism']])
		except KeyError:
			pass
	assert(len(organisms)==1)
	db.save_binomial(genid,list(organisms)[0])
	newrec = annotate(records,preserve_anno=True)
	assid = db.make_assembly(newrec,genid)
	for contig in newrec:
		db.save_from_record(assid,contig)
	return assid








def slurp_genbank(gbfilename):
	pass

def slurp_fasta(fastafilename):
	pass





if __name__=="__main__":
	pass