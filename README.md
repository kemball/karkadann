# karkadann
The Unicorn Progect 2

A complete rewrite of the GCF software and pipeline.

-Prodigal.py-

Annotate takes a genbank file and runs it through prodigal to reannotate all the genes. This is a bit rude considering genbank files mostly contain gene annotations already, but they can be patchy or missing or just inconsistent.


*annotate(gbfile)*
Returns a SeqRecord, requires a filename. 



