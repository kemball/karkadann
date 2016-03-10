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

#####`get_cursor()`

Returns a cursor the database. It's thread safe, I promise. Use like so:

>with get_cursor as cur:
>
>	#database actions

If an exception is thrown database stuff will not be committed. In fact database actions are not committed until the with block is exited. Don't be afraid to open a bunch, though, they're cheap. _they're not, but they will be_

#####`Genome(db_id,genome_name)`

Creates a new Genome object. If created with a db_id, populates itself from the database. If created with a genome_name, represents an unsaved Genome. 

* `.save()` Saves the genome to the database.
* `.delete()` Deletes the genome from the database. (Handles registered binomial names)
* `.is_real()` Checks if the Genome has been saved to the database, and if so, returns its id.
* `.added()` Returns the DateTime the genome was added to the database.
* `.binomial(binomial=None)` Acts a flexible getter/setter for the binomial name associated with this genome.

Genome.fetch(genome_name) fetches and populates a Genome object from the database and returns it.

####`Assembly(db_id,records,genome_id,)`

Creates a new Assembly object. if created with a db_id, populates itself. If created with a records iterable and a genome object, makes a new unsaved Assembly object.

* `.save()` Saves the assembly to the database.
* `.delete()` Deletes the assembly from the database. (Handles the gbfile backup)
* `.is_real()` Checks if the Assembly has been saved to the database, and if so, returns its id.
* `.is_record()` Returns an iterator over the gbfile underneath.
* `.save_record()` Helper function to generate portable filenames that aren't used yet for saving the genbank records.

#####`save_contig(assembly_id,sequence,accession=None)`

Saves contigs into the database with optional accession. Sequence is a string. `str(record.seq)` is what you want.

#####`save_from_record(assembly_id,record)`

Exactly what it sounds like. Saves the sequence as a contig and the CDS features as genes.



#####`save_genes(contig_id,features)`

Saves a bunch of genes to the specified contig, using features with locations and so on.

#####`save_gene(contig_id,translation,start,end,strand,accession=None)

Saves a single gene to the contig.

##hmm.py

This should have the hmm-parsing and database-updating code, but currently NIY.

##slurp.py

Tools for importing genbank files. Deciding on useful names, validating accession numbers, species names, etc. All that logic goes in here so it doesn't have to go anywhere else.

#####`annotate_save(gbfile,genomename)`

Reads a genbank file and a name, makes a genome for it, makes an assembly, saves the record, saves all the contigs with accessions as long as the record.id is the accession number, and annotates it. 