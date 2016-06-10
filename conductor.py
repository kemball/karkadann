import multiprocessing as mp
import subprocess as sp
from time import time
import karkadann
from karkadann.promer import promer_score
from karkadann.database import get_cursor, Cluster,Gene,Assembly,Contig,domain_max,num_cores
from karkadann.assimilate import assimilate_from_ncbi,assimilate_from_fasta
from karkadann.hmm import scan_assembly
from karkadann.cluster_call import call_clusters
from karkadann.domains import assign_groups, ortho_score
from karkadann.uclust import calc_domain_max
from itertools import combinations

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--cores',default= num_cores,type=int,help="number of cores to use. Default all.")


inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('--fasta', action="store_true",help="Specify a file is in fasta format")
inputtype.add_argument('--genbank',action="store_true",help="Specify a file is in genbank format")
parser.add_argument('--genome',nargs='+',type=str,help="A list of paths to the genomes to import. ")

parser.add_argument("--scan",action="store_true",help="scan all unscanned genomes")

parser.add_argument("--call",action="store_true",help="call gene clusters in scanned genomes")

parser.add_argument("--orthogroup",action="store_true", help="do allvall blast and declare orthogroups for genes."
                                                             " 'cluster' to just do genes in an identified cluster")
parser.add_argument("--promer",action="store_true",help="calculate all uncalculated promer scores for clusters.")

parser.add_argument("--uclust",choices = karkadann.cluster_call.types_of_clusters+["all"],help="calculate the domain_max scores"
                                                                                       " for a specified class of cluster. ")
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
		parser.error("What kind of genome file is that? I'm not equipped to guess.")
	p.close()
elif args.scan:
	with get_cursor() as cur:
		cur.execute("select distinct a.id from assemblies as a where a.id not in "
		            "(select distinct assembly_id from contigs join genes on genes.contig=contigs.id "
		            "join hits on hits.gene = genes.id);")
		assems = Assembly.get_many([x for (x,) in cur.fetchall()])
	assems = list(assems)
	p = mp.Pool(processes= max(args.cores//3,len(assems)))
	p.map(scan_assembly,assems)
	p.close()

elif args.call:
	with get_cursor() as cur:
		cur.execute("select distinct contigs.id from contigs join genes on genes.contig=contigs.id "
					" join hits on hits.gene =genes.id join clusters on clusters.gene= genes.id "
					"where clusters.id is null and hits.id is not null;")
		contigs = list(Contig.get_many([x for (x,) in cur.fetchall()]))
	print "%s uncalled contigs detected, calling clusters in them..." % len(contigs)
	p = mp.Pool(processes=args.cores)
	p.map(call_clusters,contigs)
	p.close()
elif args.orthogroup:
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
elif args.promer:
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
elif args.uclust:
	if args.uclust =="all":
		for clusttype in karkadann.cluster_call.types_of_clusters:
			calc_domain_max(clusttype)
	else:
		calc_domain_max(args.uclust)


files = sp.check_output("ls /home/kemball/diatom/plum/*.gb", shell=True).split()



def doroghazi_metric(clustera,clusterb):
	total =  ortho_score(clustera,clusterb)+2.0*domain_max(clustera,clusterb)+promer_score(clustera,clusterb)
	return total/4.0
