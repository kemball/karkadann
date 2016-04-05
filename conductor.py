import multiprocessing as mp
import subprocess as sp
from time import time


p = mp.Pool(4,maxtasksperchild=10)

files = sp.check_output("ls /home/kemball/diatom/actinobacteria_class/genbank/*.gb",shell=True).split()
files = files[10:15]

from karkadann.assimilate import assimilate_from_ncbi
before = time()
genomes = p.map(assimilate_from_ncbi,files) #genomes is a map object
print "assimilation in parallel takes %s seconds" % (time()-before)

assems = [x.assemblies().next() for x in genomes] #should support iteration tho

from karkadann.hmm import scan_assembly
before = time()
list(p.map(scan_assembly,assems))
print "parallel scanning takes %s secodns" % (time()-before)


def contig_flat(assemblylist):
	for assembly in assemblylist:
		for contig in assembly.contigs():
			yield contig

from karkadann.cluster_call import call_clusters
listcontigs = list(contig_flat(assems))
print len(listcontigs)
before = time()
clusts = p.map(call_clusters,listcontigs)
print "parallel mapping takes %s seconds" %(time()-before)
