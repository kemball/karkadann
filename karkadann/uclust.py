import os
import tempfile
import subprocess as sp
from collections import defaultdict
import itertools
from time import time


from karkadann.database import Cluster,Hit,config,domain_max,save_domain_max_many

usearch_location = config.get('usearch','usearch_location')
if not os.access(usearch_location,os.R_OK|os.X_OK):
	print "Can't find usearch executable at %s" % usearch_location
	raise IOError


def uclust_all(allhits,identity):
	if not len(allhits):
		return {}
	hitclusters = defaultdict(list)

	with tempfile.NamedTemporaryFile(suffix='.fas') as ntfasta:
		before = time()
		allhits.sort(key = lambda x:-len(x.seq))
		print "sorting takes %s " % (time()-before)
		# hitseqs should be sorted by length or clustering breaks.
		for hit in allhits:
			ntfasta.write(">%s\n%s\n" %(hit.is_real(),hit.seq))
		ntfasta.flush()
		with tempfile.NamedTemporaryFile(suffix='.uc') as tempuc:
			sp.call([usearch_location,"--cluster",ntfasta.name,"--uc",tempuc.name,"--id",str(identity),"--usersort"])
			for line in tempuc:
				if line[0] not in "SH":
					continue
				linedata = line.split('\t')
				clusternum = linedata[1]
				hitid = int(linedata[-2])

				hitclusters[clusternum].append(hitid)

	returnclusters = {}

	for key in hitclusters.keys():
		returnclusters[key] = list(Hit.get_many(hitclusters[key]))

	return returnclusters


corehmms = dict(nrps="Condensation",PKS_I="PKS_KS",PKS_II="PKS_KS",PKS_III="PKS_KS")


def calc_domain_max(cluster_kind):
	if cluster_kind not in corehmms.keys():
		raise ValueError("I don't know how to deal with cluster class %s" % cluster_kind)
	clusters = list(Cluster.by_kind(cluster_kind))
	domain_max_dict = {}
	for gca,gcb in itertools.combinations(clusters,2):
		# flip through all the clusters, if we have domain_max scores for all of them, just return that.
		# if we don't, recalculate them.
		curscore = domain_max(gca,gcb)
		if curscore is None:
			break
		else:
			domain_max_dict[(gca,gcb)] = curscore

	else:
		return domain_max_dict
	relevanthits = []
	hitsbyclust = defaultdict(list)
	for c in clusters:
		for g in c.gene_list:
			for h in g.hits():
				if h.hmm == corehmms[cluster_kind] and h.score>30:
					relevanthits.append(h)
					hitsbyclust[c].append(h)

	domain_max_dict = {(a,b):0.0 for (a,b) in itertools.combinations(clusters,2)}
	for identity in [x/10.0 for x in xrange(0,11)]:
		hc = uclust_all(relevanthits,identity)
		hit_to_uclust = {}
		for c in hc.keys():
			for h in hc[c]:
				# This needs to be indexed by db_id because comparison between
				# Hit objects is finicky and will only work if they're literally the same Hit instance.
				# I wish there was a way to mark them immutable.
				hit_to_uclust[h.is_real()]=c
		tosave = []
		for gca,gcb in domain_max_dict.keys():
			score = 0.0

			ahits = hitsbyclust[gca]
			bhits = hitsbyclust[gcb]

			if len(ahits) == 0 or len(bhits) == 0:
				continue

			aclusts = [hit_to_uclust[ah._id] for ah in ahits]
			bclusts = [hit_to_uclust[bh._id] for bh in bhits]

			for ac in aclusts:
				for bc in bclusts:
					if ac == bc:
						score += 1

			score /=(len(ahits) + len(bhits)+0.0)
			if score>=.5:
				domain_max_dict[(gca,gcb)] = identity
				tosave.append((gca,gcb,identity))
		save_domain_max_many(tosave)
	return domain_max_dict






def _set_of_domains(clustera):
	domaintypes = defaultdict(list)
	for gene in clustera.genes():
		for hit in gene.hits():
			domaintypes[hit.hmm].append(hit)
	return domaintypes


