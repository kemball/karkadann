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
		m = re.match(">?(\S+)_(\d+) # (\d+) # (\d+) # (-?\d+) # ID=([^;])+",rec.description)
		if m:
			name,id_number,start,end,strand,prod_id = m.groups()
			start = int(start)
			end = int(end)
			strand = int(strand)
			location = SeqFeature.FeatureLocation(start,end,strand)
			sequence = str(rec.seq)
			qualifiers = {'translation':sequence}
			# multiple features go on the same record. This returns the name to keep track of what goes where.
			feature = SeqFeature.SeqFeature(location=location,type="CDS",strand=strand,id=id_number,qualifiers=qualifiers)
			yield name,feature


def annotate(gbfile):
	record = list(SeqIO.parse(gbfile,'gb'))
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
	for contigrec in record:
		if contigrec.features and contigrec.features[0].type=="source":
			sf = contigrec.features[0]
			#the id here is the accession number
			contigrec.features = [sf] + features_of_contigs[contigrec.id]
		else:
			contigrec.features = features_of_contigs[contigrec.id]
	return record






if __name__=="__main__":
	betteranno = annotate('/home/kemball/actinobacteria_class/genbank/120971.gb')
	SeqIO.write(betteranno,'test.gb','genbank')
