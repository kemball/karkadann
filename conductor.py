import multiprocessing as mp
import subprocess as sp
from time import time
from karkadann.promer import promer_score
from karkadann.database import get_cursor, Cluster
from karkadann.assimilate import assimilate_from_ncbi
from karkadann.hmm import scan_assembly
from karkadann.cluster_call import call_clusters

p = mp.Pool(maxtasksperchild=10)

files = sp.check_output("ls /home/kemball/diatom/actinobacteria_class/genbank/*.gb", shell=True).split()
files = files[0:5]


before = time()
genomes = p.map(assimilate_from_ncbi, files)  # genomes is a map object
print "assimilation in parallel takes %s seconds" % (time() - before)

assems = [x.assemblies().next() for x in genomes]  # should support iteration tho


before = time()
list(p.map(scan_assembly, assems))
print "parallel scanning takes %s seconds, or %s seconds per" % ((time() - before), (time() - before) / len(assems))


def contig_flat(assemblylist):
	for assembly in assemblylist:
		for contig in assembly.contigs():
			yield contig



listcontigs = list(contig_flat(assems))
print len(listcontigs)
before = time()
clusts = p.map(call_clusters, listcontigs)
print "parallel cluster calling takes %s seconds or %s seconds per" % (
(time() - before), (time() - before) / len(listcontigs))
p.close()

def splat_promer(args):
		return promer_score(*args)

np = mp.Pool(maxtasksperchild=100)

from itertools import combinations
with get_cursor() as cur:
	cur.execute("select distinct classification from clusters;")
	clusterlist = cur.fetchall()
for (clustertype,) in clusterlist:
	with get_cursor() as cur:
		arglist = []
		cur.execute("select distinct id from clusters where classification=%s;", (clustertype,))
		clustsbytype = [Cluster(db_id=x) for (x,) in cur]
	for (ca,cb) in combinations(clustsbytype,2):
			arglist.append((ca, cb))

	print "packed up %s clusters for clustertype %s" % (len(arglist),clustertype)
	if len(arglist>10):
		np.map(splat_promer,arglist)
	else:
		map(splat_promer, arglist)
