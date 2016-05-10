from tempfile import mkdtemp,NamedTemporaryFile
import os
from shutil import rmtree
from fcntl import lockf, LOCK_EX, LOCK_UN
from Bio import SeqIO
from database import config, data_location
from database import save_orthogroup, most_recent_batch, start_batch
from database import Gene
from database import num_cores
import subprocess as sp
import re
from time import time
from collections import defaultdict

orthomcl_location = config.get('orthomcl', 'orthomcl_location')
user = config.get('orthomcl', 'username')
password = config.get('orthomcl', 'password')
db = config.get('orthomcl', 'dbname')
mcl_location = sp.check_output('which mcl', shell=True).strip()


# has to be sandboxed
def cycle_blast(input="goodProteins.fasta"):
	format_proc = sp.Popen("formatdb -t cluster -i %s -p T -n cluster"%input,shell=True)
	assert os.path.exists('cluster.pin')

	records = SeqIO.parse(input,'fasta')
	chunks = [[]]*num_cores
	for i,rec in enumerate(records):
		chunks[i%num_cores].append(rec)

	for number,chunk in enumerate(chunks):
		if len(chunk):
			SeqIO.write(chunk,'goodProteins%s.fasta'%number,'fasta')
	format_proc.wait()
	procs =[]
	for number,chunk in enumerate(chunks):
		if not len(chunk):
			continue
		np = sp.Popen("blastp -db cluster -query goodProteins%s.fasta -out orthoallvall%s.txt -outfmt 6 -num_alignments 10000 -evalue 1e-5 "%(number,number),shell=True)
		procs.append(np)
	map(lambda x: x.wait(),procs)
	sp.call("cat orthoallvall* > orthoallvall.txt",shell=True)
	return "orthoallvall.txt"


def _call_mcl(prot_records):
	homedir = os.getcwd()
	for prot in prot_records:
		prot.description = ""
		prot.name = ""
	mcl_sandbox = mkdtemp()
	try:
		os.chdir(mcl_sandbox)
		os.mkdir("fasta")
		SeqIO.write(prot_records, 'fasta/prots.fasta', 'fasta')
		sp.call("%s/bin/orthomclFilterFasta fasta 10 20" % orthomcl_location, shell=True)
		assert os.path.exists('goodProteins.fasta')
		sp.call("formatdb -t cluster -i goodProteins.fasta -p T -n cluster", shell=True)
		assert os.path.exists('cluster.pin')

		before = time()
		cycle_blast(input='goodProteins.fasta')
		print "blastp step of orthomcl takes %s seconds for %s proteins" %(time()-before,len(prot_records))
		assert os.path.exists('orthoallvall.txt')
		sp.call("%s/bin/orthomclBlastParser orthoallvall.txt fasta >> ss.txt" % orthomcl_location, shell=True)
		assert os.path.exists('ss.txt')
		before = time()
		f = open(os.path.join(data_location, 'lock.lock'), 'w')
		lockf(f, LOCK_EX)
		f.write(str(os.getpid()))
		print "waited on orthomcl lock for %s seconds" % (time() - before)
		try:
			sp.call(
				"%s/bin/orthomclInstallSchema %s/orthomcl.config orthomcl.log" % (orthomcl_location, orthomcl_location),
				shell=True)
			sp.call('%s/bin/orthomclLoadBlast %s/orthomcl.config ss.txt'
			        % (orthomcl_location, orthomcl_location), shell=True)
			sp.call('%s/bin/orthomclPairs %s/orthomcl.config orthomcl.log cleanup=yes'
			        % (orthomcl_location, orthomcl_location), shell=True)
			sp.check_output('%s/bin/orthomclDumpPairsFiles %s/orthomcl.config 1>out.log 2>error.log'
			                % (orthomcl_location, orthomcl_location), shell=True)
			assert os.path.exists('mclInput')
			sp.call('%s mclInput --abc -I 1.5 -o mclOutput 2>error.log' % mcl_location, shell=True)
			assert os.path.exists('mclOutput')
			sp.call('%s/bin/orthomclMclToGroups a 1 < mclOutput > groups.txt'
			        % (orthomcl_location), shell=True)
			assert os.path.exists('groups.txt')

		finally:
			deleteleftovertables = '''\
			drop table CoOrtholog;\
			drop table InParalog;\
			drop table Ortholog;\
			drop table SimilarSequences;\
			'''
			sp.call("mysql -e' %s' -u %s -p%s %s" % (deleteleftovertables, user, password, db), shell=True)
			lockf(f, LOCK_UN)
			f.close()
		with open('groups.txt', 'r') as returnf:
			return returnf.read()
	finally:
		rmtree(mcl_sandbox)
		os.chdir(homedir)


def parse_groups(grouptext, batch=None):
	if not batch:
		batch = most_recent_batch()
	for line in grouptext:
		m = re.match('([^:]+):', line)
		if m:
			orthogroup = m.group()
		else:
			continue
		genes = re.findall('(\d+)\|\d+_\d+', line)
		for gene in genes:
			gene_id = int(gene)
			save_orthogroup(Gene.get(gene_id), orthogroup, batch)


def assign_groups(genes):
	# just a giant iterator of gene objects, please.
	newbatch = start_batch()
	prots = [g.record for g in genes]
	if not len(prots):
		return
	gtext = _call_mcl(prots)
	parse_groups(gtext, batch=newbatch)


def ortho_score(clustera, clusterb, batch=None):
	if not batch:
		batch = most_recent_batch()
	total_genes = len(clustera.gene_list)+len( clusterb.gene_list)
	if total_genes <= 2:
		print "cluster a %s and cluster b %s have not enough genes between them. " % (clustera.is_real(),clusterb.is_real())
		return 0
	orthogroups_a = defaultdict(int)
	for gene in clustera.gene_list:
		orthogroups_a[gene.orthogroup(batch=batch)] += 1
	if "" in orthogroups_a.keys():
		raise ValueError("A gene in cluster %s has no orthogroup"% clustera.is_real())
	orthogroups_b = defaultdict(int)
	for gene in clusterb.gene_list:
		orthogroups_b[gene.orthogroup(batch=batch)] += 1
	if "" in orthogroups_b.keys():
		raise ValueError("A gene in cluster %s has no orthogroup"% clusterb.is_real())
	allkeys = list(set(orthogroups_a.keys()) | set(orthogroups_b.keys()))
	sames = 0
	for key in allkeys:
		sames += min(orthogroups_a[key], orthogroups_b[key])
	# The -2/-2 here is to correct for the shared orthogroup
	# that conceivably is always present: the cluster-starting hit.
	return (sames * 2.0-2.0) / (total_genes - 2.0)
