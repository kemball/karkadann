import os
import subprocess as sp
from threading import Thread
from tempfile import NamedTemporaryFile as ntf
from database import hmm_location, Hit
from Bio import SearchIO, SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC


def _call_hmmer(hmm, inputproteins):
	with ntf(prefix="/dev/shm/") as inputfasta:
		with ntf(prefix="/dev/shm/") as hmmoutput:
			SeqIO.write(inputproteins, inputfasta.name, 'fasta')
			hmmfile = os.path.join(hmm_location, hmm+'.hmm')
			sp.call(['hmmsearch', '-o', hmmoutput.name, hmmfile, inputfasta.name])
			QRS = SearchIO.parse(hmmoutput, format="hmmer3-text")
			for qr in QRS:
				for hit in qr:
					yield hit.id, hit.bitscore


def profile(genes, hmm):
	def aa(genes):
		for gene in genes:
			yield SeqRecord.SeqRecord(Seq.Seq(gene.translation, IUPAC.protein), id=str(gene.is_real()))

	for gene_id, score in _call_hmmer(hmm, aa(genes)):
		new_hit = Hit(gene=gene_id, score=score, hmm=hmm)
		if score > 0:
			new_hit.save()


import os


def list_hmms():
	fnames = os.listdir(hmm_location)
	return [name[:-4] for name in fnames if name.endswith('.hmm')]


def scan_assembly(assembly):
	# holy threading batman
	contigs = assembly.contigs()
	hmms = list_hmms()
	skein = []
	for c in contigs:
		for hmm in hmms:
			yarn = Thread(target=profile, args=(c.genes(), hmm), )
			yarn.start()
			skein.append(yarn)
	for yarn in skein:
		yarn.join()
