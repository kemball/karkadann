
import os
import subprocess as sp
import threading
from tempfile import NamedTemporaryFile as ntf
from database import hmm_location
from Bio import SearchIO,SeqIO,SeqRecord
from database import Gene


def call_hmmer(hmm,inputproteins):
	with ntf(prefix="/dev/shm/") as inputfasta:
		with ntf(prefix="/dev/shm/") as hmmoutput:
			SeqIO.write(inputproteins,inputfasta,'fasta')
			hmmfile = os.path.join(hmm_location, hmm)
			sp.call(['hmmsearch','-o',hmmoutput.name, hmmfile, inputfasta.name])
			hits = SearchIO.parse(hmmoutput, format="hmmer3")
			for x in hits:
				print x

if __name__=="__main__":
	input = SeqRecord.SeqRecord("QLEDKNS"*100)
	call_hmmer("AfsA.hmm",[input])