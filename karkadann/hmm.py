
import os
import subprocess as sp
import threading

from Bio import SearchIO

#call hmmer

def call_hmmer(hmmfile,inputproteins):
	sp.call(['hmmsearch',hmmfile,inputproteins])
	SearchIO.parse(hmmfile,format=hmmer3-text)
	


