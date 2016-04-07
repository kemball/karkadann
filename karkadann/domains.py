from tempfile import mkdtemp
import os
from shutil import rmtree
from fcntl import lockf,LOCK_EX,LOCK_UN
from Bio import SeqIO
from database import config,data_location
from database import save_orthogroup,most_recent_batch,start_batch
import subprocess as sp
import re
from time import time

orthomcl_location = config.get('orthomcl','orthomcl_location')
user  = config.get('orthomcl','username')
password = config.get('orthomcl','password')
db = config.get('orthomcl','dbname')
mcl_location = sp.check_output('which mcl',shell=True).strip()


def _call_mcl(prot_records):
	homedir = os.getcwd()
	for prot in prot_records:
		prot.description = ""
		prot.name = ""
	mcl_sandbox = mkdtemp()
	try:
		os.chdir(mcl_sandbox)
		os.mkdir("fasta")
		SeqIO.write(prot_records,'fasta/prots.fasta','fasta')
		sp.call("%s/bin/orthomclFilterFasta fasta 10 20" % orthomcl_location, shell=True)
		assert os.path.exists('goodProteins.fasta')
		sp.call("formatdb -t cluster -i goodProteins.fasta -p T -n cluster",shell=True)
		assert os.path.exists('cluster.pin')
		# need to thread this or hack up the input into pieces or something.
		# what is this 2005
		sp.call('blastp -db cluster -query goodProteins.fasta -out orthoallvall.txt -outfmt 6 -num_alignments 10000 -evalue 1e-5 -num_threads 8',shell=True)
		assert os.path.exists('orthoallvall.txt')
		sp.call("%s/bin/orthomclBlastParser orthoallvall.txt fasta >> ss.txt"%orthomcl_location,shell=True)
		assert os.path.exists('ss.txt')
		before = time()
		f = open(os.path.join(data_location,'lock.lock'), 'w')
		lockf(f, LOCK_EX)
		f.write(str(os.getpid()))
		print "waited on orthomcl lock for %s seconds" % (time()-before)
		try:
			sp.call("%s/bin/orthomclInstallSchema %s/orthomcl.config orthomcl.log" %(orthomcl_location,orthomcl_location),shell=True)
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
			deleteleftovertables='''\
			drop table CoOrtholog;\
			drop table InParalog;\
			drop table Ortholog;\
			drop table SimilarSequences;\
			'''
			sp.call("mysql -e' %s' -u %s -p%s %s"%(deleteleftovertables,user,password,db),shell=True)
			lockf(f, LOCK_UN)
			f.close()
		with open('groups.txt','r') as returnf:
			return returnf.read()
	finally:
		rmtree(mcl_sandbox)
		os.chdir(homedir)

def parse_groups(grouptext,batch=None):
	if not batch:
		batch = most_recent_batch()
	for line in grouptext:
		m = re.match('([^:]+):',line)
		if m:
			orthogroup = m.group()
		else:
			continue
		genes = re.findall('(\d+)\|\1',line)
		for gene in genes:
			gene_id = int(gene)
			save_orthogroup(gene_id,orthogroup,batch)

def assign_groups(genes):
	# just a giant iterator of gene objects, please.
	newbatch = start_batch()
	prots = [g.record for g in genes]
	gtext = _call_mcl(prots)
	parse_groups(gtext,batch=newbatch)
