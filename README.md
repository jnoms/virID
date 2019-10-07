# virID
Viral Identification and Discovery - A viral characterization pipeline built in [Nextflow](https://www.nextflow.io/).

## Description
The purpose of virID is to assemble and classify in next-generation sequencing data. The steps of this pipeline are as follows:
1) Input fastqs are split into a paired file containing interleaved paired reads, and an unpaired file containing unpaired reads.
2) Reads are assembled with the SPAdes assembler - SPAdes, metaSPAdes, or rnaSPAdes can be specified.
3) bwa mem is used to map reads back to contigs.
4) The DIAMOND and megablast aligners are used for taxonomic assignment of assembled contigs.
5) DIAMOND and megablast output files are translated to a taxonomic output following last-common-ancestor calculation for each query contig.
6) (virID v2.0+) Contigs are queried with megablast against a nonredundant database of common cloning vectors. Contigs that are assigned to these sequences are flagged.
6) DIAMOND and megablast taxonomic outputs and contig count information are merged to a comprehensive taxonomic output, and unassigned contigs are flagged. Counts outputs, one for each DIAMOND and megablast, are generated which display the counts at every taxonomic level.

While this pipeline is organized in Nextflow, every process is capabale of being used and tested independently of Nextflow as a bash or python script. Execute any python script (with -h) or bash script (without arguments) for detailed instructions on how to run them.

## Prepare programs and databases
1. Download Nextflow version 19.07.0. I recommend using the [Conda package manager](https://docs.conda.io/en/latest/miniconda.html).
`conda install nextflow=19.07.0`

2. Download dependencies. I recommend generating a conda environment using virID_environment.yml, available in this github.
`conda env create -f virID_environment.yml`. This environment will contain all tools and packages required for running this pipeline. You can test that everything works independently by directly calling each shell and python script.

3. Configure the databases:
This workflow requires three databases: 1) A DIAMOND database for search, 2) A nucleotide blast database for search, and 3) a contaminant nucleotide blast database (virID v2.0+).

A DIAMOND database can be generated by following the DIAMOND [manual](https://github.com/bbuchfink/diamond/raw/master/diamond_manual.pdf). I recommend generating a database using the RefSeq nonredundant protein fastas present at [ftp://ftp.ncbi.nlm.nih.gov/refseq/release/].

For blast seach, you must download the blast v5 nucleotide database using the same version of blast+ used in this script and installed via the conda .yml file. `update_blastdb.pl --blastdb_version 5 nr_v5 --decompress`. There is sometimes a bug where this database will not work - you may have to delete the last line of nt_v5.nal that specifies a volume of the blast database that is not actually present.

A nucleotide contaminant blast database is included with this pipeline. This database was generated from [Univec_ core](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/?), which is a nonredundant list of common laboratory cloning vectors. I have also added some additional vectors to this database.


## Setting inputs and parameters
In general, all input values and parameters for this script must be entered in the nextflow.config file.  

#### Input and output
**params.out_dir:** Desired output directory.  
**params.reads:** A glob detailing the location of reads to be processed through this pipeline. NOTE, input reads for one sample should be in one fastq. This script will process the fastq into an interleaved (paired) fastq and a separate unpaired fastq. For example, "raw_reads/*fastq" is appropriate. Make sure to wrap in quotes!  

#### SPAdes
**params.spades_type:** Options are 'meta', 'rna', and 'regular'. This specifies whether to use metaspades.py, rnaspades.py, or spades.py.  
**params.temp_dir:** This specifies the location of the SPAdes temporary directory.  
**params.spades_min_length:** The minimum contig length that is output from the SPAdes process. Contigs below this length will be filtered out.  

#### DIAMOND
**params.diamond_database:** Path to the diamond database. Wrap all paths in quotes!  
**params.diamond_evalue:** The maximum evalue of a match to be reported as an alignment by DIAMOND. In general I am a fan of setting this at 10, which is quite high. Lower values, such as 0.001, are more stringent and have fewer false positives.  
**params.diamond_outfmt:** This dictates the output format from DIAMOND. This pipeline requires an outfmt of "6 qseqid stitle sseqid staxids evalue bitscore pident length", which is specified in the config file.  

#### BLAST
**params.blast_database:** Path to the blast nt database. Wrap all paths in quotes!  
**params.blast_evalue:** The maximum evalue of a match to be reported as an alignment by blast. In general I am a fan of setting this at 10, which is quite high. Lower values, such as 0.001, are more stringent and have fewer false positives.  
**params.blast_type:** This controls the 'task' switch of blast, and specified whether to run 'megablast', 'dc-megablast', or 'blastn'. I recommend 'megablast', else it will be quite slow.  
**params.blast_outfmt:** This dictates the output format from BLAST. This pipeline requires an outfmt of "6 qseqid stitle sseqid staxid evalue bitscore pident length", which is specified in the config file.  
**params.blast_log_file:** This is a file that will contain log information from the blast bash script.  
**params.blast_max_hsphs:** This is the maximum number of alignments per query-subject pair. Recommend setting at 1.  
**params.blast_max_targets:** This is the maximum number of separate hits that are allowed for each query.  
**params.blast_restrict_to_taxids:** This parameter lets you limited the blast search to a specified taxonID. This causes an extreme speedup, and is helpful for testing this pipeline. Not compatible with params.blast_ignore_taxids.  
**params.blast_ignore_taxids:** This parameter lets you ignore all hits of a particular taxonID.  


#### BLAST/DIAMOND conversion
**params.within_percent_of_top_score:** When finding the LCA of all matches for a given query sequence, this details how close to the maximum bitscore a match must be to be considered in the LCA classification. If this is set at 1, for example, all potential alignments within 1 percent of the highest bitscore for a query sequence will be considered in the LCA classification. **NOTE**: This is limited intrinsically by the DIAMOND -top parameter, which is set at 1. Thus, DIAMOND will only output assignments within 1% of the top bitscore anyway. I will add a switch to change the DIAMOND -top parameter in a future release.  
**params.taxid_blacklist:** Path to a file containing taxonIDs to be blacklisted. I have included a file in this github repository. Assignments containing one of these taxonIDs will be discarded before LCA calculation.  
**params.diamond_readable_colnames:** These are the more-readable column names that will be reported in the output from DIAMOND. If you change the outfmt, change this line accordingly.  
**params.blast_readable_colnames:** These are the more-readable column names that will be reported in the output from BLAST. If you change the outfmt, change this line accordingly.  
**params.taxonomy_column:** This details which of the colnames contains the taxonID.  
**params.score_column:** This details which of the colnames should be used for calculating the top score. I use bitscore, but you could technically set this as pident or length or evalue to sort by one of those parameters instead.  

#### BLAST/DIAMOND conversion
**params.blast_contaminant_database:** Path the the vector contaminant database.  
**params.blast_contaminant_evalue:** Required evalue for assignment. I recommend setting this to be fairly stringent.  
**params.blast_contaminant_outfmt:** Outformat. For ease, I default to the normal blast outfmt.  
**params.blast_contaminant_log_file:** Path to the log file.  
**params.blast_contaminant_type:** The -task blast parameter. I recommend megablast or else it will be way too slow, and unnecessary.  
**params.blast_contaminant_max_hsphs:** This is the maximum number of alignments per query-subject pair. Recommend setting at 1.  
**params.blast_contaminant_max_targets:** This is the maximum number of separate hits that are allowed for each query. Because this is a screen, and we only need to know if a query has at least one assignment in the contaminant database, I recommend setting at 1.  

## Configure executor and resources
The Nextflow executor, explained [here](https://www.nextflow.io/docs/latest/executor.html), dictates how Nextflow will run each process. virID is currently set up to use a SLURM cluster, but you can easily change this by altering the executor in nextflow.config. Nextflow takes care of all cluster submissions and automatically parallelizes everything.

I have set up virID with *dynamic resources!* Essentially, each process will request a different amount of resources depending on the size of the input files. In addition, upon failure the process will restart with twice the resources.

## Description of outfiles
*in progress*
