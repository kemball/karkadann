import os
import tempfile
import subprocess as sp
from collections import defaultdict
import itertools
import re
from time import time

from Bio import SeqIO

from karkadann.database import Cluster,Gene,Hit,config

usearch_location = config.get('usearch','usearch_location')
if not os.access(usearch_location,os.R_OK|os.X_OK):
	print "Can't find usearch executable at %s" % usearch_location
	raise IOError


def score_uclust_clusters(clustera, clusterb):
	a_by_domain = _set_of_domains(clustera)
	b_by_domain = _set_of_domains(clusterb)
	for hmm in a_by_domain.keys():
		ahits = a_by_domain[hmm]
		bhits = b_by_domain[hmm]
		with tempfile.NamedTemporaryFile() as ntf:
			pass


# needs to take only the hits.
def uclust_all(genes,identity):
	# need to keep track of source cluster
	allhits = []
	hitclusters = defaultdict(list)
	for gene in genes:
		for hit in gene.hits():
			allhits.append(hit)
	with tempfile.NamedTemporaryFile(suffix='.fas') as ntfasta:
		before = time()
		allhits.sort(key = lambda x:len(x._seq))
		print "sorting takes %s " % (time()-before)
		# hitseqs should be sorted by length or clustering breaks.
		for hit in allhits:
			ntfasta.write(">%s\n%s\n" %(hit._id,hit._seq))
		with tempfile.NamedTemporaryFile(suffix='.uc') as tempuc:
			sp.call([usearch_location,"--cluster",ntfasta.name,"--uc",tempuc.name,"--id",str(identity),"--usersort"])
			for line in tempuc:
				if line[0] =="#":
					continue
				linedata = line.split('\t')
				clusternum = linedata[1]
				hitid = linedata[-2]
				hitclusters[clusternum].append(hitid)
	return hitclusters


def doroghazi(list_of_clusters):
	return  # a dictionary of scores defined pairwise


def _set_of_domains(clustera):
	domaintypes = defaultdict(list)
	for gene in clustera.genes():
		for hit in gene.hits():
			domaintypes[hit.hmm].append(hit)
	return domaintypes


