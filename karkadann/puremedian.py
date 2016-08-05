
from Bio.SubsMat import MatrixInfo
from Bio import pairwise2
import numpy as np
from munkres import Munkres,print_matrix
from karkadann.database import Cluster


_gapscores = MatrixInfo.blosum62

def seq_identity(seqa,seqb):
	aligns = pairwise2.align.globaldx(seqa,seqb,_gapscores)
	a,b,score,begin,end = aligns[0]
	matches = sum([ac == bc for (ac,bc) in zip(a,b)])
	return matches / (1.0*len(a))

key_hmms = {"PKS_KS", "bt1fas", "ft1fas", "hg1D", "hg1E", "fabH", "tra_KS", "PKS_AT", "t2ks", "ene_KS", "mod_KS",
            "hyb_KS", "itr_KS", "t2fas", "t2clf","Chal_sti_synt_C","Chal_sti_synt_N","Condensation",
            "AMP_binding", "A_OX","Terpene_synth_C","phytoene_synt","Lycopene_cycl","terpene_cyclase","NapT7","fung_ggpps",
            "fung_gpps2","dmat","trichodiene_synth","LANC_like","Land_dehyd_N","Lant_dehyd_C","TIGR03731","BLS",
            "CAS","tabtoxin","strH_like","strK_like1","strK_like2","neoL_like","DOIS","valA_like","spcFG_like",
            "spcDK_like_cou","salQ","IucA_IucC","ectoine_synt","AfsA","indsynth","LipM","LipU","LipV","ToyB","TunD",
            "pur6","pur10","nikJ","nikO","MoeO5","melC","ppm","nos","lan_c_jd","terpene_jd","TOMM"
            }

def pure_median_identity(clustera,clusterb):
	ahits = {g:g.hits() for g in clustera.gene_list}
	bhits = {g:g.hits() for g in clusterb.gene_list}
	aseqs = []
	bseqs = []
	for gene,hitlist in zip(ahits.keys(),ahits.values()):
		aseqs.extend([hit.seq for hit in hitlist if hit.hmm in key_hmms])
	for gene, hitlist in zip(bhits.keys(),bhits.values()):
		bseqs.extend([hit.seq for hit in hitlist if hit.hmm in key_hmms])
	if not len(aseqs) or not(len(bseqs)):
		return 0.0
	prematrix = np.zeros((len(aseqs),len(bseqs)))
	for ai,aseq in enumerate(aseqs):
		for bi,bseq in enumerate(bseqs):
			si = seq_identity(aseq,bseq)
			prematrix[ai][bi] = si
			print si
	npm = np.array(prematrix)
	print npm
	m = Munkres()
	indices = m.compute(npm)
	print indices
	identities = [prematrix[ai][bi] for (ai,bi) in indices]
	print identities
	return np.median(identities)
