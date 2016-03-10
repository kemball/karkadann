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

Contains all of the database handling things, trying to keep all the database code in one place. I'm considering refactoring to an OO approach but not for now. Also defines/checks the data directory from config 

#####`get_cursor()`

Returns a cursor the database. It's thread safe, I promise. Use like so:

>with get_cursor as cur:
>
>	#database actions

If an exception is thrown database stuff will not be committed. In fact database actions are not committed until the with block is exited. Don't be afraid to open a bunch, though, they're cheap.

#####`save_record(record,salt=None)`

Saves the record provided in genbank format to the data_location directory. Returns the gb_record string used to look it up in the filesystem, so you can copy the data directory and update the config and everything will work correctly. If you provide a salt it'll try and use it, but if it's already in use it'll pick a new one. 

#####`read_record(assem_id)`

Returns the SeqRecord associated with that assembly id. 

#####`make_genome(name)`

Create a new genome named `name`. 

#####`get_genome(name)`

Returns the genome id associated with `name`.

#####`make_assembly(record,genome_id,kwargs)`

Creates an assembly linked to a genome. Automagically saves the record to the filesystem in the data directory.

#####`save_contig(assembly_id,sequence,accession=None)`

Saves contigs into the database with optional accession. Sequence is a string. `str(record.seq)` is what you want.

#####`save_from_record(assembly_id,record)

Exactly what it sounds like. Saves the sequence as a contig and the CDS features as genes.

#####`save_binomial(genome_id,genus_species)`

Saves a binomial name into the database for that genome.

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