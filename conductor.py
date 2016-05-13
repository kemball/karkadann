import multiprocessing as mp
import subprocess as sp
from time import time
from karkadann.promer import promer_score
from karkadann.database import get_cursor, Cluster,Gene,Assembly,domain_max
from karkadann.assimilate import assimilate_from_ncbi,assimilate_from_fasta
from karkadann.hmm import scan_assembly
from karkadann.cluster_call import call_clusters
from karkadann.domains import assign_groups, ortho_score
from itertools import combinations

p = mp.Pool(maxtasksperchild=10)

files = sp.check_output("ls /home/kemball/diatom/plum/*.gb", shell=True).split()


before = time()
genomes = p.map(assimilate_from_ncbi, files)  # genomes is a map object
print "assimilation in parallel takes %s seconds" % (time() - before)

before = time()
assems = []
for genome in genomes:
	assems.append(next(genome.assemblies()))
print "making list of assemblies took %s" % (time() - before)


before = time()
tinypool = mp.Pool(processes=2)
# scan_assembly is internally threaded.
list(tinypool.map(scan_assembly, assems))
print "parallel scanning takes %s seconds, or %s seconds per" % ((time() - before),
                                                                 (time() - before) / len(assems))
# this only does cluster-based genes which is a bit cheatz


def contig_flat(assemblylist):
	for assembly in assemblylist:
		for contig in assembly.contigs():
			yield contig

listcontigs = contig_flat(assems)

before = time()

clusts = p.map(call_clusters, listcontigs)
print "parallel cluster calling takes %s seconds " % (time() - before)


with get_cursor() as cur:
	cur.execute("select distinct gene from clusters;")
	allgenesinvolved = Gene.get_many([x for (x,) in cur.fetchall()])
before = time()
assign_groups(allgenesinvolved)
print "assigning groups takes %s seconds" % (time() - before)

p.close()


def splat_promer(args):
	x = promer_score(*args)


def splat_domain(args):
	x = ortho_score(*args)

np = mp.Pool(maxtasksperchild=10)

#
from uclust import corehmms,calc_domain_max
for clustertype in ["nrps","PKS_I","PKS_II","PKS_III"]:
	calc_domain_max(clustertype)
	with get_cursor() as cur:
		arglist = []
		cur.execute("select distinct id from clusters where classification=%s;", (clustertype,))
		results = cur.fetchall()
	clustsbytype = list(Cluster.get_many([x for (x,) in results]))
	if len(clustsbytype) < 2:
		continue
	# well we can't very well compare a cluster to itself now can we.
	for (ca, cb) in combinations(clustsbytype, 2):
		arglist.append((ca, cb))

	print "packed up %s cluster pairs for clustertype %s" % (len(arglist), clustertype)

	np.map(splat_promer, arglist)
	np.map(splat_domain, arglist)


def doroghazi_metric(clustera,clusterb):

	total =  ortho_score(clustera,clusterb)+2.0*domain_max(clustera,clusterb)+promer_score(clustera,clusterb)
	return total/4.0
