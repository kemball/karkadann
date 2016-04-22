from Bio import SeqIO
from Bio import SeqFeature,BiopythonWarning
from Bio.Alphabet import IUPAC
import warnings
from database import Assembly, Genome, Contig, Gene
from prodigal import annotate
import re
import os

from random import sample
from string import ascii_lowercase

list_of_culture_collections = [
	"DSM",
	"NRRL",
	"ATCC",
	"SCGC"
]


def _slug(text, aggressiveness=2):
	# check for allcaps+numbers words?
	# those are usually strain names...
	# Really should put culture collection hints into the names table. *shrug*.
	# TODO put culture collection hitns into the names table.
	if aggressiveness > 7:
		return _slug(text, aggressiveness=7) + "".join(sample(ascii_lowercase, 5))
	m = re.search(r'\s([A-Z0-9]+)\s', text)
	if m and aggressiveness == 2 and m.group().strip() not in list_of_culture_collections:
		return m.group().strip()
	words = text.split()
	return "".join([word[:aggressiveness] for word in words[:min(aggressiveness,4)]])




def make_standard(records):
	# ok, I know this is weird, but I'm sick of runs of Ns not being labeled and partial proteins
	# having CDS but no 'translation'. 
	for record in records:
		for feat in record.features:
			if feat.type == "CDS" and "translation" not in feat.qualifiers.keys():
				n_seq = feat.extract(record.seq)
				while (len(n_seq) % 3) > 0:
					n_seq += 'N'
				# there's a way to fetch the translation table from the record
				# but I'd hate to rely on that. Another day.
				feat.qualifiers["translation"] = [n_seq.translate(table=11)]
		m = re.search(r'NN+', str(record.seq))
		if m:
			gap_location = SeqFeature.FeatureLocation(start=m.start(), end=m.end(), strand=0)
			gap_feat = SeqFeature.SeqFeature(location=gap_location, type='assembly_gap')
			record.features.append(gap_feat)


from _mysql_exceptions import IntegrityError


def _save_unique(description):
	agg=2
	genome_name = _slug(description, agg)
	newgenome = Genome(genome_name=genome_name)
	while not newgenome.is_real():
		try:
			newgenome.save()
		except IntegrityError:
			# This is fairly spammy. I should change it when I add
			# automagical culture collection detection.
			print "failed to save genome with name %s trying something else" % newgenome._name
			agg += 1
			newgenome = Genome(genome_name=_slug(description, agg))
	return newgenome


def _save_contigs(reannotated_record,newassem):
	for record in reannotated_record:
		newcontig = Contig(seq=str(record.seq), assembly=newassem, accession=record.id)
		newcontig.save()
		gene_iterable = []
		for feat in record.features:
			if feat.type == "CDS":
				newgene = Gene(translation=feat.qualifiers['translation'][0],
				               contig=newcontig,
				               start=int(feat.location.start),
				               end=int(feat.location.end),
				               strand=int(feat.location.strand),
				               accession=feat.qualifiers.get("protein_id", [None])[0])
				gene_iterable.append(newgene)
		Gene._save_many(gene_iterable)


def assimilate_from_ncbi(ncbifile):
	# This is very scary. I understand. 100% of ncbi genbank files
	# have length DBLink fields. "assembly:bleh biosample:bleh etc"
	# If the whole field is longer than 80 characters, it warns us.
	# Seems like a lot of spam to ensure you never have to linewrap
	# when displaying a genbank file, but what do I know.

	warnings.simplefilter('ignore',BiopythonWarning)

	ncbirecord = list(SeqIO.parse(ncbifile, 'genbank'))
	make_standard(ncbirecord)
	reannotated_record = annotate(ncbirecord, preserve_anno=True)
	desc = reannotated_record[0].description
	newgenome = _save_unique(desc)

	# this better be right...
	newgenome.binomial(reannotated_record[0].annotations['organism'])
	try:
		m = re.findall(r'Assembly:(\S+)', reannotated_record[0].dbxrefs[0])
		if m:
			assem_acc = m[0]
		else:
			# they don't always have an assembly accession number.
			assem_acc = None

		newassem = Assembly(record=reannotated_record, genome=newgenome, accession=assem_acc)
		# I'd wrap this is a try/except but what could I even do? If this fails Here There Be Problems
		newassem.save()
		print "assembly for genome %s saved" % newgenome._name
		_save_contigs(reannotated_record,newassem)
		print "all contigs for genome %s saved" %newgenome._name

	except:
		newgenome.delete()
		raise
	return newgenome

#TODO
# Need a doroghazi importer


def assimilate_from_fasta(fastafile):
	# raw contig importer
	ncbirecord = list(SeqIO.parse(fastafile,'fasta',alphabet=IUPAC.IUPACAmbiguousDNA()))
	make_standard(ncbirecord)
	ng = _save_unique(os.path.basename(fastafile))
	try:
		reannotated_record = annotate(ncbirecord)
		newassem = Assembly(record=reannotated_record, genome=ng)
		newassem.save()
		_save_contigs(reannotated_record,newassem)
	except:
		ng.delete()
		raise
	return ng


if __name__ == "__main__":
	from karkadann.database import data_location
	import os

	assimilate_from_ncbi(os.path.join(data_location, 'test/testassem.gb'))
