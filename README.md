# virID
Viral Identification and Discovery - A viral characterization pipeline built in Nextflow.

## Description
The purpose of virID is to assemble and classify in next-generation sequencing data. The steps of this pipeline are as follows:
1) Input fastqs are split into a paired file containing interleaved paired reads, and an unpaired file containing unpaired reads.
2) Reads are assembled with the SPAdes assembler - SPAdes, metaSPAdes, or rnaSPAdes can be specified.
3) bwa mem is used to map reads back to contigs.
4) The DIAMOND and megablast aligners are used for taxonomic assignment of assembled contigs.
5) DIAMOND and megablast output files are translated to a taxonomic output following last-common-ancestor calculation for each query contig.
6) DIAMOND and megablast taxonomic outputs and contig count information are merged to a comprehensive taxonomic output, and unassigned contigs are flagged. Counts outputs, one for each DIAMOND and megablast, are generated which display the counts at every taxonomic level.

While this pipeline is organized in Nextflow, every process is capabale of being used and tested independently of Nextflow as a bash or python script. Run any python script (with -h) or bash script (without arguments) for detailed instructions on how to run them.

## Usage Instructions
*in progress*
