from Bio import SeqIO
from Bio import SeqFeature
from database import Assembly,Genome,Contig,Gene 
from prodigal import annotate
import re

from random import sample
from string import ascii_lowercase
def _slug(text,aggressiveness=2):
	#check for allcaps+numbers words?
	#those are usually strain names...
	if aggressiveness>5:
		return _slug(text,aggressiveness=3)+sample(ascii_lowercase,5)
	m = re.search(r'\s([A-Z0-9]+)\s',text)
	if m and aggressiveness == 1:
		return m.group().strip()
	words = text.split()
	return "".join([word[:aggressiveness] for word in words[:aggressiveness]])

import re
def make_standard(records):
	# ok, I know this is weird, but I'm sick of runs of Ns not being labeled and partial proteins
	# having CDS but no 'translation'. 
	for record in records:
		for feat in record.features:
			if feat.type=="CDS" and "translation" not in feat.qualifiers.keys():
				n_seq = feat.extract(record.seq)
				while (len(n_seq) %3)>0:
					n_seq += 'N'
				#there's a way to fetch the translation table from the record
				#but I'd hate to rely on that. Another day.
				feat.qualifiers["translation"] = [n_seq.translate(table=11)]
		m = re.search(r'NN+',str(record.seq))
		if m:
			gap_location = SeqFeature.FeatureLocation(start=m.start(),end=m.end(),strand=0)
			gap_feat = SeqFeature.SeqFeature(location=gap_location,type='assembly_gap')
			record.features.append(gap_feat)


from _mysql_exceptions import IntegrityError
def assimilate_from_ncbi(ncbifile):

	ncbirecord= list(SeqIO.parse(ncbifile,'genbank'))
	#can thread here, but rule #1 of optimization
	make_standard(ncbirecord)
	reannotated_record = annotate(ncbirecord)
	agg = 1
	desc = reannotated_record[0].description
	genome_name = _slug(desc,agg)
	newgenome = Genome(genome_name=genome_name)
	while not newgenome.is_real():
		try:
			newgenome.save()
		except IntegrityError:
			agg += 1
			newgenome = Genome(genome_name=_slug(desc,agg))

	#this better be right...
	newgenome.binomial(reannotated_record[0].annotations['organism'])

	try:
		m = re.findall(r'Assembly:(\S)+',reannotated_record[0].dbxrefs[0])
		if m:
			assem_acc = m[0]
		else:
			raise Exception("Assembly dbxref not found in ncbi-generated genbank.")

		newassem = Assembly(record=reannotated_record,genome=newgenome,accession=assem_acc)
		#I'd wrap this is a try/except but what could I even do? If this fails Here There Be Problems
		newassem.save()
		for record in reannotated_record:
			newcontig = Contig(seq=str(record.seq),assembly=newassem,accession=record.id)
			newcontig.save()
			for feat in record.features:
				if feat.type =="CDS":
					newgene = Gene(translation=feat.qualifiers['translation'],
									contig = newcontig,
									start=feat.location.start,
									end=feat.location.end,
									strand=feat.location.strand,
									accession = feat.qualifiers.get("protein_id",[None])[0]
									)
					newgene.save()
	finally:
		newgenome.delete()



if __name__=="__main__":
	from karkadann.database import data_location
	import os
	assimilate_from_ncbi(os.path.join(data_location,'test/testassem.gb'))

	
