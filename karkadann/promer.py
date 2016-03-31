import subprocess as sp
from tempfile import mkdtemp
from Bio import SeqIO
import os
import re
from shutil import rmtree


def _parse_delta(deltafilename):
	for line in open(deltafilename, 'r'):
		if re.match(">(\S+) (\S+) (\d+) (\d+)", line):
			leftname, rightname, leftlen, rightlen = re.match(">(\S+) (\S+) (\d+) (\d+)", line).groups()
			leftlen, rightlen = int(leftlen), int(rightlen)
			leftAA = [0] * leftlen
			rightAA = [0] * rightlen
		m = re.match("(\d+) (\d+) (\d+) (\d+) \d+ \d+ \d+", line)
		if m:
			leftstart, leftend, rightstart, rightend = m.groups()
			leftstart = int(leftstart)
			leftend = int(leftend)
			rightstart = int(rightstart)
			rightend = int(rightend)
			# probably requires explanation. Assigning into tuples is great
			leftstart, leftend = min(leftstart, leftend), max(leftstart, leftend)
			rightstart, rightend = min(rightstart, rightend), max(rightstart, rightend)
			leftAA[leftstart:leftend] = [1] * (leftend - leftstart)
			rightAA[rightstart:rightend] = [1] * (rightend - rightstart)
	return sum(leftAA) / (2.0 * leftlen) + sum(rightAA) / (2.0 * leftlen)


def promer_score(clustera, clusterb):
	return _call_promer(clustera.fna(), clusterb.fna())


def _call_promer(fasta1, fasta2):
	promer_sandbox = mkdtemp()
	try:
		os.chdir(promer_sandbox)
		SeqIO.write(fasta1, 'one.fna', 'fasta')
		SeqIO.write(fasta2, 'two.fna', 'fasta')
		sp.call(["promer --mum --prefix=toast one.fna two.fna 2>mummerout.txt"], shell=True)
		return _parse_delta("toast.delta")
	finally:
		rmtree(promer_sandbox)
