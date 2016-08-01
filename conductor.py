import multiprocessing as mp
import subprocess as sp
from time import time
import karkadann
from karkadann.promer import promer_score
from karkadann.database import *
from karkadann.assimilate import assimilate_from_ncbi, assimilate_from_fasta, assimilate_from_antismash
from karkadann.hmm import scan_assembly
from karkadann.cluster_call import call_clusters
from karkadann.domains import assign_groups, ortho_score
from karkadann.uclust import calc_domain_max
from itertools import combinations

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--cores', default=num_cores, type=int, help="number of cores to use. Default all.")
parser.add_argument('--force', action="store_true", help="Try extra hard to do things.")

inputtype = parser.add_mutually_exclusive_group()
inputtype.add_argument('--fasta', action="store_true", help="Specify a file is in fasta format")
inputtype.add_argument('--genbank', action="store_true", help="Specify a file is in genbank format")
parser.add_argument('--genome', nargs='+', type=str, help="A list of paths to the genomes to import. ")

parser.add_argument("--antismash", nargs='+', type=str, help="A directory containing antismash results to import.")

parser.add_argument("--scan", action="store_true", help="scan all unscanned genomes")

parser.add_argument("--call", action="store_true", help="call gene clusters in scanned genomes")

parser.add_argument("--orthogroup", choices=["cluster", "all"],
                    help="do allvall blast and declare orthogroups for genes."
                         " 'cluster' to just do genes in an identified cluster")
parser.add_argument("--promer", action="store_true", help="calculate all uncalculated promer scores for clusters.")

parser.add_argument("--uclust", choices=karkadann.cluster_call.types_of_clusters + ["all"],
                    help="calculate the domain_max scores"
                         " for a specified class of cluster. ")

parser.add_argument("--table", choices=["clusters", "genomes"],
                    help="prints a tab-separated table of information about what you asked for.")

parser.add_argument("--network", choices=["D", "all"],
                    help="prints a tsv file suitable for cytoscape import of the metric specified")

parser.add_argument("--export", nargs='+', type=str, help="Print a cluster in the specified format to stdout")

parser.add_argument("--type", choices=karkadann.cluster_call.types_of_clusters,
                    help="filter results by a specific cluster type. Applies to table and network output")

args = parser.parse_args()

if args.cores:
	karkadann.database.num_cores = args.cores

with get_cursor() as cur:
	cur.execute("set @max_connections = 10000;")

if args.genome:
	p = mp.Pool(processes=args.cores)
	if args.fasta:
		genomes = p.map(assimilate_from_fasta, args.genome)
	elif args.genbank:
		genomes = p.map(assimilate_from_ncbi, args.genome)
	else:
		parser.error("What kind of genome file is that? I'm not prepared to guess.")
	p.close()
if args.antismash:
	for adir in args.antismash:
		assimilate_from_antismash(adir)
if args.scan:
	before = time()
	with get_cursor() as cur:
		if not args.force:
			cur.execute("select distinct a.id from assemblies as a where a.id not in "
			            "(select distinct assembly_id from contigs join genes on genes.contig=contigs.id "
			            "join hits on hits.gene = genes.id);")
		else:
			print "Re-scanning all genomes"
			cur.execute("delete from hits;")
			cur.execute("select distinct a.id from assemblies;")
		assems = Assembly.get_many([x for (x,) in cur.fetchall()])
	assems = list(assems)
	p = mp.Pool(processes=min(args.cores // 3, len(assems)+1) )
	p.map(scan_assembly, assems)
	p.close()
	print "%s assemblies scanned in %s seconds." % (len(assems), before - time())

if args.call:
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
		print "%s contigs detected, re-calling clusters in them..." % len(contigs)
	p = mp.Pool(processes=args.cores)
	p.map(call_clusters, contigs)
	p.close()
	print "Completed successfully in %s seconds." % (time() - before)
if args.orthogroup:
	before = time()
	if args.orthogroup == "cluster":
		with get_cursor() as cur:
			cur.execute("select distinct gene from clusters;")
			genes = Gene.get_many([x for (x,) in cur.fetchall()])
	elif args.orthogroup == "all":
		with get_cursor() as cur:
			cur.execute("select id from genes;")
			genes = Gene.get_many([x for (x,) in cur.fetchall()])
	else:
		genes = []
		parser.exit("Which genes should I calculate orthogroups for? options are 'cluster' and 'all'")
	assign_groups(genes)
	print "assigned orthogroups in %s seconds" % (before - time())
if args.promer:
	before = time()
	if not args.type:
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters;")
			clusters = list(Cluster.get_many([x for (x,) in cur.fetchall()]))
	else:
		with get_cursor() as cur:
			cur.execute("select distinct id from clusters where classification=%s;", (args.type,))
			clusters = list(Cluster.get_many([x for (x,) in cur.fetchall()]))


	def splat_promer(args):
		return promer_score(*args)


	p = mp.Pool(processes=args.cores)
	arglist = []
	for ca, cb in combinations(clusters, 2):
		arglist.append((ca, cb))
	p.map(splat_promer, arglist)
	p.close()
	print "%s promer scores calculated in %s seconds" % (len(arglist), time() - before)
if args.uclust:
	before = time()
	if args.uclust == "all":
		for clusttype in karkadann.cluster_call.types_of_clusters:
			typebefore = time()
			calc_domain_max(clusttype)
			print "cluster type %s took %s seconds." % (clusttype, time() - typebefore)
	else:
		typebefore = time()
		calc_domain_max(args.uclust)
		print "cluster type %s took %s seconds." % (args.uclust, time() - typebefore)
	print "Finished uclustering successfully. %s seconds." % (time() - before)
if args.table:
	if args.table == "clusters":
		if not args.type:
			with get_cursor() as cur:
				cur.execute("select distinct id from clusters;")
				clusters = Cluster.get_many([x for (x,) in cur.fetchall()])
		else:
			with get_cursor() as cur:
				cur.execute("select distinct id from clusters where classification=%s", (args.type,))
				clusters = Cluster.get_many([x for (x,) in cur.fetchall()])
		print "identity\tname\tkind"
		for clust in clusters:
			print "\t".join([str(clust._id), clust.name, clust._kind])
	elif args.table == "genomes":
		with get_cursor() as cur:
			cur.execute("select distinct id from genomes;")
			genomes = Genome.get_many([x for (x,) in cur.fetchall()])
		print "identity\tname\tgenus species"
		for g in genomes:
			print "\t".join([str(g._id), g._name, g.binomial()])
if args.network:

	def doroghazi_metric(clustera, clusterb):
		total = ortho_score(clustera, clusterb) + 2.0 * domain_max(clustera, clusterb) + promer_score(clustera,
		                                                                                              clusterb)
		return total / 4.0


	if args.network == "D":
		if not args.type:
			parser.exit("A type is required to output D metrics.")
		clusters = Cluster.by_kind(args.type)
		print "\t".join(
			["caid", "caname", "metric", "cbid", "cbname", "D", "orthoscore", "domain_maxscore", "promerscore"])


		def threadlet((ca, cb)):
			oscore = ortho_score(ca, cb)
			dmaxscore = domain_max(ca, cb)
			pscore = promer_score(ca, cb)
			row = [ca._id, ca.name, "D-metric", cb._id, cb.name, doroghazi_metric(ca, cb), oscore, dmaxscore, pscore]
			# if oscore>=.5 and dmaxscore>=.7 and pscore >=.5:
			return "\t".join(map(str, row))


		p = mp.Pool(processes=args.cores)
		for result in p.imap_unordered(threadlet, combinations(clusters, 2), chunksize=100):
			if result:
				print result
if args.export:
	if args.type and args.type == "fasta":
		with get_cursor() as cur:
			names = ",".join(args.export)
			cur.execute("select id from cluster_names where name in (%s);", (names,))
			clusts = Cluster.get_many([x for (x,) in cur.fetchall()])
		for clust in clusts:
			print clust.fna()
	else:
		# genbank format
		with get_cursor() as cur:

			query = "select id from cluster_names where name in (%s);" % (",".join(["%s" for x in args.export]))
			cur.execute(query, args.export)
			clusts = Cluster.get_many([x for (x,) in cur.fetchall()])
		clusts = list(clusts)
		for clust in clusts:
			print clust.genbank()

def families(group):
	allclusts = list(Cluster.by_kind(group))
	families = [[x] for x in allclusts]
	def cutoffpass(clustera,clusterb):
		oscore = ortho_score(clustera,clusterb)
		if oscore < .5:
			return False
		pscore = promer_score(clustera,clusterb)
		if pscore <.5:
			return False
		dscore = domain_max(clustera,clusterb)
		if dscore <.7:
			return False
		return True
	shouldstop = False
	while not shouldstop:
		shouldstop = True
		singletons = [x[0] for x in families if len(x)==1]
		for sin in singletons:
			for family in families:
				for member in family:
					if cutoffpass(member,sin):
						families.remove([sin])
						family.append(sin)
						shouldstop = False
	return families


