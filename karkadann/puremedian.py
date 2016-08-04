
from Bio.SubsMat import MatrixInfo
from Bio import pairwise2
from karkadann.database import Cluster


_gapscores = MatrixInfo.blosum62

def seq_identity(seqa,seqb):
	aligns = pairwise2.align.globaldx(seqa,seqb,_gapscores)
	a,b,score,begin,end = aligns[0]
	matches = sum([ac == bc for (ac,bc) in zip(a,b)])
	return matches / (1.0*len(a))

def pure_median_identity(clustera,clusterb):
	ahits = {g:g.hits() for g in clustera.gene_list}
	bhits = {g:g.hits() for g in clusterb.gene_list}
