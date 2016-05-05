import os
import re
import subprocess as sp
from tempfile import NamedTemporaryFile as ntf
from threading import Thread,Semaphore

from Bio import SearchIO, SeqIO, SeqRecord, Seq
from Bio.Alphabet import IUPAC

from database import hmm_location, Hit, Gene, num_cores


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
				# qr.sort()
				# I'm kind of hoping this sorts by hit strength.
				# worth checking. I guess it doesn't matter anyway.

				for hit in qr:
					scores[hit.id] = max(scores[hit.id], hit.bitscore)
					for hsp in hit.hsps:
						def appropriate_hyphens(m):
							return '-' * len(m.group(0))

						if len(hsp.hit.seq) > 100:
							hitseq = re.sub('PPPPP+', appropriate_hyphens, str(hsp.hit.seq))
							hitseq = hitseq.translate(None,'-*').upper()
							yield hit.id, hsp.bitscore, hitseq


def profile(genes, hmms):
	# TODO(MAYBE) return meaningful results from this
	def aa(geneset):
		for gene in geneset:
			yield SeqRecord.SeqRecord(Seq.Seq(gene.translation, IUPAC.protein), id=str(gene.is_real()))

	protset = list(aa(list(genes)))
	skein = []

	def threadlet(miniprot,hmm,semaphore):
		semaphore.acquire()
		hit_list = []
		for gene_id, score,sequence in _call_hmmer(hmm, miniprot):
			new_hit = Hit(gene=Gene.get(gene_id), score=score,seq=sequence, hmm=hmm)
			if score > 0:
				hit_list.append(new_hit)
		Hit._save_many(hit_list)
		semaphore.release()
	sem = Semaphore(num_cores)
	for hmm in hmms:
		yarn = Thread(target=threadlet,args=(protset,hmm,sem))
		yarn.start()
		skein.append(yarn)
	for yarn in skein:
		yarn.join()



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
