# karkadann
The Unicorn Progect 2

A complete rewrite of the GCF software and pipeline.

##Install

`python setup.py install` should install karkadann to your current python environment. I recommend virtenvs and pip for dependencies and so on, but MySQLdb is legendarily difficult to install. It's full of pitfalls, so you're always going to have some work to do unless you have it already. Anyway.

`python setup.py test` will run automatic tests of all the code that exists so far. Some of it relies on testassem.gb being in a particular place but that's being phased out.

##Usage

I've documented all the function signatures here but the API is still half-baked and subject to change. At some point there'll be a set of functions that are available for export and guaranteed not to change, but that point is not today.

###Prodigal.py

Annotate takes a SeqIO record iterator and runs it through prodigal to reannotate all the genes. This is a bit rude considering genbank files mostly contain gene annotations already, but they can be patchy or missing or just inconsistent.


#####`annotate(gbrecord,preserve_anno=False)`
Returns a prodigated gbrecord. If preserve_anno is True, it'll only add prodigal genes that don't overlap with genbank annotations. 


###Database.py

Contains all of the database handling things, also defines/checks the data directory from config.

####`get_cursor()`

Returns a cursor the database. It's thread safe, I promise. Use like so:

>with get_cursor as cur:
>
>	#database actions

If an exception is thrown database stuff will not be committed. In fact database actions are not committed until the with block is exited. 





####`Genome(genome_name)`

Creates a new Genome object. If created with a genome_name, represents an unsaved Genome. 

* `.get(db_id)` Fetches a Genome object from the database by id.
* `.save()` Saves the genome to the database.
* `.delete()` Deletes the genome from the database. (Handles registered binomial names)
* `.is_real()` Checks if the Genome has been saved to the database, and if so, returns its id.
* `.added()` Returns the DateTime the genome was added to the database.
* `.binomial(binomial=None)` Acts a flexible getter/setter for the binomial name associated with this genome.

`Genome.fetch(genome_name)` fetches and populates a Genome object from the database and returns it.

#### `Assembly(records,genome_id,)`

Creates a new Assembly object.  If created with a records iterable and a genome object, makes a new unsaved Assembly object.

* `.get(db_id)` Fetches an Assembly object from the database by id.
* `.save()` Saves the assembly to the database.
* `.delete()` Deletes the assembly from the database. (Handles the gbfile backup)
* `.is_real()` Checks if the Assembly has been saved to the database, and if so, returns its id.
* `.is_record()` Returns an iterator over the gbfile underneath.
* `.save_record()` Helper function to generate portable filenames that aren't used yet for saving the genbank records. You shouldn't need this at all.
* `.contigs()` Returns a generator of contigs in this assembly.

`Assembly.fetch(genome_id)` fetches a list of assemblies based on the genome_id. 99% of the time there won't be more than one but the database supports that on purpose, so it doesn't make sense to ignore it here.

####`Contig(seq,assembly,accession)`

* `.get(db_id)` Fetches a Contig from the database by id.
* `.save()` Saves the Contig
* `.delete()` Deletes the contig from the database. (handles genes correctly)
* `.is_real()` Checks if the Contig has been saved to the database, and if so, returns its id.
* `.seq()` Returns the sequence of the contig.
* `.acc()` Returns the accession number for the contig.
* `.genes()` Returns a generator of the Genes in the contig.

####`Gene(translation,contig,start,end,strand,accession)`

* `.get(db_id)` Fetches a Gene from the database by id.
* `.save()` Saves the Gene to the database.
* `.delete()` Deletes the gene from the database.
* `.is_real()` Checks if the Gene has been saved, and if so returns its id.
* `.location` Returns the FeatureLocation for the Gene.
* `.translation` Returns the AA sequence for the Gene.
* `.orthogroup` Returns a string representing the orthogroup this gene belongs to. A property. `.orthogroup(batch=batchnum)` will return the orthogroup string for that particular batch.

####`Hit(Gene,score,hmm)`

* `.get(db_id)` Fetches a Hit from the database by id.
* `.save()` Saves the hit to the database.
* `.delete()` Deletes the hit from the database.
* `.is_real()` Checks if the hit has been saved, and if so returns its id.
* `.score` Pseudoproperty. Returns the bitscore of the hit.

####`Cluster(gene_list,classification)`
Reads a cluster from the database or makes a new one from a classification and a gene_list. The id isn't numeric but a string based on its classification and the number of existing clusters.

* `.get(db_id)` Fetches a cluster object from the database
* `.save()` Saves the cluster to the database.
* `.delete()` Deletes the cluster from the database.
* `.is_real()` Checks if the cluster has been saved, and if so returns its id.
* `.fna()` Returns a SeqRecord object representing the underlying DNA for the cluster
* `.faa()` Returns a list of SeqRecords for each protein in the cluster in order.


#### Batch functions 
* `most_recent_batch()` Returns the id for the most recent orthomcl batch created. The default, usually.
* `start_batch()` Creates a new batch for orthomcl processing.
* `save_orthogroup(gene,orthogroup,batch=most_recent_batch())` Labels a particular gene as belonging to a particular orthogroup for a particular batch.

You probably shouldn't need these, they're called internally by domains.assign_groups, but you never know.


##hmm.py
#### `list_hmms()`
 Lists all the hmms available in the data directory karkadann maintains.
#### `scan_assembly(assembly)`
Runs all the hmms against the specified assembly. Totally unthreaded(for now).
#### `profile(genes,hmm)`
Takes an iterator of Gene objects and an hmm string 'AfsA.hmm' for example. Does the needful(scans the hmms against the gene translations and then saves the results into the database).

##cluster_call.py
#### `_classify(gene)`
Returns the string classifying a gene into a seed for a gene cluster. You almost certainly shouldn't need to call this.

#### `call_clusters(contig)`
Calls all the gene clusters in a given contig and puts them into the database.

##promer.py
#### `promer_score(clusterone,clustertwo)`
Returns the average involvement of the two gene clusters. Tidies up after itself, promer is messy. Parallelizes great, though. Checks if the promer score has already been calculated. If it has, fetches that value from the database. If it hasn't, calculates it and stores it in the database.

##domains.py
#### `assign_groups`(genes)
Automatically and safely creates a new 'orthomcl batch'. All orthomcl group callings are referenced by batch number so that a second orthomcl clustering can be done without clobbering the first. 

#### `domain_score`(clustera,clusterb,batch)
Returns the shared-domain score for two clusters.

##assimilate.py

Tools for importing genbank files. Deciding on useful names, validating accession numbers, species names, etc. All that logic goes in here so it doesn't have to go anywhere else.

####`assimilate_from_ncbi(ncbifile)`

Does exactly what it sounds like. Takes just one minute for a 10Mb genome, counting time taken for prodigal, merging the features together, writing all the stuff to the database. Prodigal takes 20-40 seconds no matter what. Picks a genome name from the description of the genbank file. Reads the assembly accession number off DBLINK. Could cross-reference against Biosample in the future if that is required or interesting.

####`assimilate_from_fasta(fastafile)`

Does exactly what it sounds like. Takes a fasta file of contigs and tries to make educated guesses about what the organism name etc is likely to be. Mostly gets it wrong, but it tries. Annotates and saves to the database.

Will write another helper that slurps data out of the gcf database.

####`make_standard(records)`

This adds 'translation' qualifiers to partial genes that hang off contig borders and detects and annotates assembly gaps. It doesn't make any real decisions but it means I don't have to worry about how it works anywhere else in the code.
