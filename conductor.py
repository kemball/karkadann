import multiprocessing as mp
import subprocess as sp
from time import time
import karkadann
from karkadann.promer import promer_score
from karkadann.database import *
from karkadann.assimilate import assimilate_from_ncbi,assimilate_from_fasta
from karkadann.hmm import scan_assembly
from karkadann.cluster_call import call_clusters
from karkadann.domains import assign_groups, ortho_score
from karkadann.uclust import calc_domain_max
from itertools import combinations

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--cores',default= num_cores,type=int,help="number of cores to use. Default all.")
parser.add_argument('--force',action="store_true",help="Try extra hard to do things.")

inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('--fasta', action="store_true",help="Specify a file is in fasta format")
inputtype.add_argument('--genbank',action="store_true",help="Specify a file is in genbank format")
parser.add_argument('--genome',nargs='+',type=str,help="A list of paths to the genomes to import. ")

parser.add_argument("--scan",action="store_true",help="scan all unscanned genomes")

parser.add_argument("--call",action="store_true",help="call gene clusters in scanned genomes")

parser.add_argument("--orthogroup",choices = ["cluster","all"], help="do allvall blast and declare orthogroups for genes."
                                                             " 'cluster' to just do genes in an identified cluster")
parser.add_argument("--promer",action="store_true",help="calculate all uncalculated promer scores for clusters.")

parser.add_argument("--uclust",choices = karkadann.cluster_call.types_of_clusters+["all"],help="calculate the domain_max scores"
                                                                                       " for a specified class of cluster. ")

parser.add_argument("--table",choices=["clusters","genomes"],help="prints a tab-separated table of information about what you asked for.")

args = parser.parse_args()

if args.cores:
	karkadann.database.num_cores = args.cores

if args.genome:
	p = mp.Pool(processes=args.cores)
	if args.fasta:
		genomes = p.map(assimilate_from_fasta,args.genome)
	elif args.genbank:
		genomes = p.map(assimilate_from_ncbi,args.genome)
	else:
		parser.error("What kind of genome file is that? I'm not prepared to guess.")
	p.close()
elif args.scan:
	before = time()
	with get_cursor() as cur:
		cur.execute("select distinct a.id from assemblies as a where a.id not in "
		            "(select distinct assembly_id from contigs join genes on genes.contig=contigs.id "
		            "join hits on hits.gene = genes.id);")
		assems = Assembly.get_many([x for (x,) in cur.fetchall()])
	assems = list(assems)
	p = mp.Pool(processes= max(args.cores//3,len(assems)))
	p.map(scan_assembly,assems)
	p.close()
	print "%s assemblies scanned in %s seconds." %(len(assems),before-time())

elif args.call:
	before = time()
	if not args.force:
		with get_cursor() as cur:
			cur.execute("select distinct contigs.id from contigs join genes on genes.contig=contigs.id "
						" join hits on hits.gene =genes.id join clusters on clusters.gene= genes.id "
						"where clusters.id is null and hits.id is not null;")
			contigs = list(Contig.get_many([x for (x,) in cur.fetchall()]))
		print "%s uncalled contigs detected, calling clusters in them..." % len(contigs)
	else:
		with get_cursor() as cur:
			cur.execute("select distinct contigs.id from contigs;")
			contigs = list(Contig.get_many([x for (x,) in cur.fetchall()]))
		for contig in contigs:
			for clust in contig.clusters():
				clust.delete()
		print "%s contigs detected, re-calling clusters in them..."
	p = mp.Pool(processes=args.cores)
	p.map(call_clusters,contigs)
	p.close()
	print "Completed successfully in %s seconds." %(time()-before)
elif args.orthogroup:
	before= time()
	if args.orthogroup=="cluster":
		with get_cursor() as cur:
			cur.execute("select distinct gene from clusters;")
			genes = Gene.get_many([x for (x,) in cur.fetchall()])
	elif args.orthogroup=="all":
		with get_cursor() as cur:
			cur.execute("select id from genes;")
			genes = Gene.get_many([x for (x,) in cur.fetchall()])
	else:
		genes = []
		parser.exit("Which genes should I calculate orthogroups for? options are 'cluster' and 'all'")
	assign_groups(genes)
	print "assigned orthogroups in %s seconds" %(before-time())
elif args.promer:
	before = time()
	with get_cursor() as cur:
		cur.execute("select distinct id from clusters;")
		clusters = list(Cluster.get_many([x for (x,) in cur.fetchall()]))

	def splat_promer(args):
		return promer_score(*args)
	p = mp.Pool(processes=args.cores)
	arglist = []
	for ca,cb in combinations(clusters,2):
		arglist.append((ca,cb))
	p.map(splat_promer,arglist)
	p.close()
	print "%s promer scores calculated in %s seconds" %(len(arglist),time()-before)
elif args.uclust:
	before = time()
	if args.uclust =="all":
		for clusttype in karkadann.cluster_call.types_of_clusters:
			typebefore = time()
			calc_domain_max(clusttype)
			print "cluster type %s took %s seconds." %(clusttype,time()-typebefore)
	else:
		typebefore = time()
		calc_domain_max(args.uclust)
		print "cluster type %s took %s seconds." %(args.uclust,time()-typebefore)
	print "Finished uclustering successfully. %s seconds."%(time()-before)
elif args.table:
	if args.table=="clusters":
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters;")
			clusters = Cluster.get_many([x for (x,) in cur.fetchall()])
		print "identity\tname\tkind"
		for clust in clusters:
			print "\t".join([str(clust._id),clust.name,clust._kind])
	elif args.table=="genomes":
		with get_cursor() as cur:
			cur.execute("select distinct id from genomes;")
			genomes = Genome.get_many([x for (x,) in cur.fetchall()])
		print "identity\tname\tgenus species"
		for g in genomes:
			print "\t".join([str(g._id),g._name,g.binomial()])
else:
	print "Sorry, I'm not sure what you're asking me to do."





def doroghazi_metric(clustera,clusterb):
	total =  ortho_score(clustera,clusterb)+2.0*domain_max(clustera,clusterb)+promer_score(clustera,clusterb)
	return total/4.0
