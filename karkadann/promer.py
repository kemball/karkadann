import subprocess as sp
from _mysql_exceptions import IntegrityError
from tempfile import mkdtemp
from Bio import SeqIO
import os
import re
from shutil import rmtree


def _parse_delta(deltafilename):
	leftAA = []
	rightAA = []
	leftlen = 1
	rightlen = 1
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
	return sum(leftAA) / (2.0 * leftlen) + sum(rightAA) / (2.0 * rightlen)


def promer_score(clustera, clusterb):
	from database import get_cursor
	ida = clustera._id
	idb = clusterb._id
	if ida==idb:
		return 1
	if not ida or not idb:
		print "raising valuerror"
		raise ValueError("scoring unreal clusters is impossible")
	# enforce ida<idb
	ida,idb = min(ida,idb),max(ida,idb)
	with get_cursor() as cur:
		cur.execute("select score from promer where l= %s and r =%s;",(ida,idb))
		result = cur.fetchone()
		if result:
			return result[0]
	with get_cursor() as cur:
		score = _call_promer(clustera.fna(),clusterb.fna())
		try:
			cur.execute("insert into promer (score,l,r) values(%s,%s,%s);",(score,ida,idb))
		except IntegrityError:
			print "possible threading problem caught. Do clusters %s and %s have a score?" %(ida,idb)
	return score



def _call_promer(fasta1, fasta2):
	promer_sandbox = mkdtemp(prefix='/dev/shm/')
	try:
		os.chdir(promer_sandbox)
		SeqIO.write(fasta1, 'one.fna', 'fasta')
		SeqIO.write(fasta2, 'two.fna', 'fasta')
		sp.call(["promer --mum --prefix=toast one.fna two.fna 2>mummerout.txt"], shell=True)
		return _parse_delta("toast.delta")
	finally:
		rmtree(promer_sandbox)
