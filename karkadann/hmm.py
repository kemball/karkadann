
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


