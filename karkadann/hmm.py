import os
import subprocess as sp
from threading import Thread
from tempfile import NamedTemporaryFile as ntf
from database import hmm_location, Hit, Gene
from Bio import SearchIO, SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC


def _call_hmmer(hmm, inputproteins):
	inputproteins = list(inputproteins)
	scores = {}
	for ip in inputproteins:
		scores[ip.id] = 0

	with ntf(prefix="/dev/shm/") as inputfasta:
		with ntf(prefix="/dev/shm/") as hmmoutput:
			SeqIO.write(inputproteins, inputfasta.name, 'fasta')
			hmmfile = os.path.join(hmm_location, hmm + '.hmm')
			sp.call(['hmmsearch', '-o', hmmoutput.name, hmmfile, inputfasta.name])
			hmmoutput.flush()
			hmmoutput.seek(0)
			QRS = SearchIO.parse(hmmoutput, format="hmmer3-text")
			for qr in QRS:
				# there's *always* a QR, even though it's usually empty.
				qr.sort()
				# I'm kind of hoping this sorts by hit strength.
				# worth checking.

				for hit in qr:
					scores[hit.id] = max(scores[hit.id], hit.bitscore)
	for ip in inputproteins:
		if scores[ip.id]:
			yield ip.id, scores[ip.id]


def profile(genes, hmms):
	def aa(geneset):
		for gene in geneset:
			yield SeqRecord.SeqRecord(Seq.Seq(gene.translation, IUPAC.protein), id=str(gene.is_real()))

	protset = list(aa(list(genes)))
	skein = []

	def threadlet(hmm,protset):
		hit_list = []
		for gene_id, score in _call_hmmer(hmm, protset):
			new_hit = Hit(gene=Gene(db_id=int(gene_id)), score=score, hmm=hmm)
			if score > 0:
				hit_list.append(new_hit)
		Hit.save_many(hit_list)
	for hmm in hmms:
		yarn = Thread(target=threadlet,args=(hmm,protset))
		yarn.start()
		skein.append(yarn)
	for yarn in skein:
		yarn.join()


import os


def list_hmms():
	fnames = os.listdir(hmm_location)
	return [name[:-4] for name in fnames if name.endswith('.hmm')]


def scan_assembly(assembly):
	# holy threading batman
	contigs = assembly.contigs()
	hmms = list_hmms()
	from itertools import chain
	all_genes = chain(*[c.genes() for c in contigs])
	profile(all_genes, hmms)
