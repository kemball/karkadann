import subprocess as sp 
from Bio import SeqIO,SeqFeature
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import re
from tempfile import NamedTemporaryFile as ntf
import os
from collections import defaultdict

#requires prodigal, but not sure how to check for that.
#i'll figure it out later

#I really need to figure out exactly the conversions done


def call_prodigal(fastafile):
	"""Invokes prodigal on a provided fasta file, returns the SeqRecord produced by -a.
	Everything is done in temporary files kept on virtual filesystem."""
	#check if file exists blah blah
	with ntf(prefix='/dev/shm/',delete=True,suffix='.prot') as protfile,ntf(prefix='/dev/shm/',delete=True,suffix='.out') as prod:
		sp.call(['prodigal','-i',fastafile,'-a',protfile.name,'-o',prod.name])
		#you can't close over temporary files, so the .parse generator can't generate once this returns
		#hence list. sucks to be you memory
		return list(SeqIO.parse(protfile.name,'fasta'))
		

# >NZ_KK070020.1_7 # 4671 # 5855 # 1 # ID=3_7;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.734
def parse_prodigal(prodigal_record):
	"""Takes the pseudo-fasta and returns a list of SeqFeatures."""
	for rec in prodigal_record:
		#each one of these records is a feature
		m = re.match(">?(\S+)_(\d+) # (\d+) # (\d+) # (-?\d+) # ID=([^;]+);",rec.description)
		if m:
			name,id_number,start,end,strand,prod_id = m.groups()
			start = int(start)
			end = int(end)
			strand = int(strand)
			location = SeqFeature.FeatureLocation(start,end,strand)
			sequence = str(rec.seq)
			qualifiers = {'translation':sequence,'prodigal_id':prod_id}
			# multiple features go on the same record. This returns the name to keep track of what goes where.
			feature = SeqFeature.SeqFeature(location=location,type="CDS",strand=strand,id=id_number,qualifiers=qualifiers)
			yield name,feature

def overlap(one,two):
	#this is designed to operate on SeqFeature locations
	#sharing a single base doesn't count. I DECREE
	if one.start == two.end or two.start == one.end:
		return False
	elif one.start in two or one.end in two:
		return True
	elif two.start in one or two.end in one:
		return True
	else:
		return False







def annotate(record,preserve_anno = False):
	record = list(record)
	#could parallelize but won't help much
	with ntf(prefix='/dev/shm/',delete=True,suffix=".fna") as fastafile:
		# prodigal can handle the genbank files
		# but in order to generate pseudofastas with accessions instead of species names
		# it has to use fasta.
		# to avoid excessive disk sadness, /dev/shm files are kept on RAM
		SeqIO.write(record,fastafile.name,'fasta')
		gene_calls = parse_prodigal(call_prodigal(fastafile.name))

	features_of_contigs = defaultdict(list)
	for name,feature in gene_calls:
		features_of_contigs[name].append(feature)
	if not preserve_anno:
		#these are entirely independent, and could be threaded. 
		#But what's the point? parallelizing calls to annotate would be faster and require
		# less effort
		for contigrec in record:
			if contigrec.features and contigrec.features[0].type=="source":
				sf = contigrec.features[0]
				#the id here is the accession number
				contigrec.features = [sf] + features_of_contigs[contigrec.id]
			else:
				contigrec.features = features_of_contigs[contigrec.id]
	else:
		for contigrec in record:
			newfeats = features_of_contigs[contigrec.id]
			oldfeats = contigrec.features
			keepfeats = []
			for n in newfeats:
				for o in oldfeats:
					#source types and assembly gaps don't count
					if o.type in ["gene","CDS"] and overlap(n.location,o.location):
						break
				else:
					#we should keep this new feature, it doesn't overlap
					keepfeats.append(n)
			contigrec.features.extend(keepfeats)

	return record






if __name__=="__main__":
	overlap_test()
	#this is nto a real test case but it seems to work so...
	testrec=SeqIO.parse('/home/kemball/bunicorn/data/test/testassem.gb','genbank')
	betteranno = annotate(testrec)
	SeqIO.write(betteranno,'test.gb','genbank')
	testrec=SeqIO.parse('/home/kemball/bunicorn/data/test/testassem.gb','genbank')
	preserved_anno = annotate(testrec,preserve_anno=True)
	SeqIO.write(preserved_anno,'test2.gb','genbank')
