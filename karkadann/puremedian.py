
from Bio.SubsMat import MatrixInfo
from Bio import pairwise2,SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np
from munkres import Munkres,print_matrix
from shutil import rmtree

import os
from tempfile import mkdtemp
import subprocess as sp


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
	with get_cursor() as cur:
		cur.execute("select score from pure_median where l = %s and r= %s;",(clustera._id,clusterb._id))
		for res in cur:
			return res[0]
		pmid = _pure_median_identity(clustera,clusterb)
		cur.execute("replace into pure_median (score,l,r) values (%s,%s,%s);",(pmid,clustera._id,clusterb._id))
		return pmid


def _pure_median_identity(clustera,clusterb):
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
	# ok, so the following looks crazy but I promise it makes quote sense unquote
	# for some reason the munkres implementation does not like matrices taller than
	# they are long. Longer than they are tall, sure, but not taller than they are long.
	# for my purposes I guess it doesn't matter, just reverse the two.
	if len(aseqs)>len(bseqs):
		aseqs,bseqs = bseqs,aseqs
	oldwd = os.getcwd()
	sandbox = mkdtemp(prefix='/dev/shm/')
	os.chdir(sandbox)
	try:
		arecords = [SeqRecord(seq,id=str(num)) for num,seq in enumerate(aseqs)]
		brecords = [SeqRecord(seq,id=str(num)) for num,seq in enumerate(bseqs)]
		SeqIO.write(arecords,'ahits.fasta',"fasta")
		SeqIO.write(brecords,'bhits.fasta',"fasta")
		bp = sp.Popen("blastp -query ahits.fasta -subject bhits.fasta -out results.txt -outfmt 6 ",shell=True)
		bp.wait()

		prematrix = np.ones((len(aseqs),len(bseqs)))
		for line in open("results.txt",'r'):
			linelist = line.split()
			qid,sid,pid = linelist[0:3]
			qid = int(qid)
			sid = int(sid)
			pid = float(pid)
			# the qid is always from a, the sid from b
			# Munkres operates on a cost matrix, so inverse identity here.
			# I guess that would be called the percent difference.
			# That's also why it's an array of ones. 100% different.
			prematrix[int(qid)][int(sid)]=1-float(pid)/100.0
	finally:
		os.chdir(oldwd)
		rmtree(sandbox)
	npm = prematrix
	m = Munkres()
	try:
		indices = m.compute(npm)
	except ValueError:
		print "computing munkres threw"
		print
		print npm
	# convert from percent difference back to identities
	identities = [1-prematrix[ai][bi] for (ai,bi) in indices]
	return np.median(identities)
