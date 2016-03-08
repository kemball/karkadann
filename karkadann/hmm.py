
import os
import subprocess as sp
import threading
from tempfile import NamedTemporaryFile as ntf
from database import hmm_location
from Bio import SearchIO,SeqIO


def call_hmmer(hmm,inputproteins):
	hmmfile=os.path.join(hmm_location,hmm)
	sp.call(['hmmsearch',hmmfile,inputproteins])
	hits = SearchIO.parse(hmmfile,format=hmmer3-text)
	for x in hits:
		print x

if __name__=="__main__":
	from database import import_hmm
	from database import data_location
	testgb  = os.path.join(data_location,"test/testassem.gb")
	from slurp import annotate_save
	assid = annotate_save(testgb,'testinghmms')
	from database import get_cursor as gc
	with gc() as cur:
		cur.execute("select genes.id,genes.translation from contigs \
			join assemblies\
				on contigs.assembly_id = assemblies.id \
			join genes on contigs.id=genes.contig where assemblies.id = %s;",(assid,))
		
	sname = import_hmm('/home/kemball/bunicorn/hmms/abmotifs.hmm')
	call_hmmer(sname,)


