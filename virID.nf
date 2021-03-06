#!/usr/bin/env nextflow

nextflow.preview.dsl=2

// Log file must be set here to access laucnhDir information
params.log_file = "${workflow.launchDir}/${params.out_dir}/reports/virID.log"

// Require nextflow version 19.10 or higher
if( !nextflow.version.matches('19.10+') ) {
    println "This workflow requires Nextflow version 19.10 or greater -- You \
    are running version $nextflow.version"
    exit 1
}

//============================================================================//
// Set up modules
//============================================================================//
include './bin/modules/process_read_pairs' params(params)

include './bin/modules/spades_assembly' params(params)

include './bin/modules/diamond'  params(params)

include blast_to_LCA as convert_diamond from './bin/modules/blast_to_LCA' params(
  LCA_top_percent: params.LCA_top_percent,
  out_dir: params.out_dir,
  column_names: params.diamond_readable_colnames,
  source: 'diamond',
  taxid_blacklist: params.taxid_blacklist,
  within_percent_of_top_score: params.within_percent_of_top_score = "1",
  taxonomy_column: params.taxonomy_column,
  score_column:params.score_column
  )

include './bin/modules/blast' params(params)

include blast_to_LCA as convert_blast from './bin/modules/blast_to_LCA' params(
  LCA_top_percent: params.LCA_top_percent,
  out_dir: params.out_dir,
  column_names: params.blast_readable_colnames,
  source: 'blast',
  taxid_blacklist: params.taxid_blacklist
  )

include blast as blast_contaminant from './bin/modules/blast' params(
  blast_database: params.blast_contaminant_database,
  blast_evalue: params.blast_contaminant_evalue,
  blast_outfmt: params.blast_contaminant_outfmt,
  blast_log_file: params.blast_contaminant_log_file,
  blast_type: params.blast_contaminant_type,
  blast_max_hsphs: params.blast_contaminant_max_hsphs,
  blast_max_targets: params.blast_contaminant_max_targets,
  out_dir: "$params.out_dir/contaminant",
  log_file: params.log_file
  )

include './bin/modules/bwa_mem' params(params)

include './bin/modules/generate_output' params(params)

include './bin/modules/fastq_to_fasta' params(params)

include get_counts as get_counts_blast from './bin/modules/get_counts' params(
  out_dir: params.out_dir,
  source: "blast"
  )

include get_counts as get_counts_diamond from './bin/modules/get_counts' params(
  out_dir: params.out_dir,
  source: "diamond"
  )

//============================================================================//
// Defining functions
//============================================================================//
def sampleID_set_from_infile(input) {
  // The purpose of this function is to take an input list of file paths and
  // to return a channel of sets of structure basename:file_path.
  sample_set = []
  for (String item : file(input)) {
    file = file(item)
    name = file.baseName
    sample_set.add([name, file])
  }
  ch = Channel.from(tuple(sample_set))
  return ch
}

//============================================================================//
// Define workflows
//============================================================================//
workflow assembly_pipeline {

  get: input_ch
  main:
  // Generate contigs
  process_read_pairs(input_ch) \
    | spades_assembly

  contigs = spades_assembly.out
    .filter{ it[1].size() > 0 }

  // Map reads back to assembly
  contigs
    .join(process_read_pairs.out) \
    | bwa_mem_contigs

  // DIAMOND processing
  diamond(contigs)
    .filter{ it[1].size() > 0 } \
    | convert_diamond

  // BLAST Processing
  blast(contigs)
    .filter{ it[1].size()>0 } \
    | convert_blast

  // Assignment against contaminant database
  blast_contaminant(contigs)

  // Merge channels and generate output
  convert_blast.out
    .join(convert_diamond.out)
    .join(contigs)
    .join(bwa_mem_contigs.out)
    .join(blast_contaminant.out) \
    | generate_output
}

workflow read_pipeline {
  get: input_ch
  main:

  // Process read pairs
  process_read_pairs(input_ch)

  // Convert input to fasta
  fastq_to_fasta(process_read_pairs.out)

  // Run DIAMOND
  if (params.reads_pipeline_no_diamond == "F") {
    diamond(fastq_to_fasta.out)
      .filter{ it[1].size() > 0 } \
      | convert_diamond \
      | get_counts_diamond
  }

  // Run BLAST
  if (params.reads_pipeline_no_blast == "F") {
    blast(fastq_to_fasta.out)
      .filter{ it[1].size()>0 } \
      | convert_blast \
      | get_counts_blast
  }

  // Run contaminant blast
  blast_contaminant(fastq_to_fasta.out)
}

//============================================================================//
// Validate inputs
//============================================================================//

if( (params.assembly_pipeline == "F") && (params.reads_pipeline == "F") ) {
  error "Either params.assembly_pipeline or params.reads_pipeline must \
  be set to T."
}
else if( (params.assembly_pipeline == "T") && (params.reads_pipeline == "T") ) {
  error "Only one of params.assembly_pipeline or params.reads_pipeline must \
  be set to T."
}

if( (params.reads_pipeline_no_diamond == "T") && (params.reads_pipeline_no_blast == "T") ) {
  error "Only one of params.reads_pipeline_no_diamond or \
  params.reads_pipeline_no_blast may be set to 'T'."
}

//============================================================================//
// Define main workflow
//============================================================================//
workflow {

  main:
    input_ch = sampleID_set_from_infile(params.reads)

    if ( params.assembly_pipeline == "T" ) {
      assembly_pipeline(input_ch)
    }
    else if ( params.reads_pipeline == "T" ) {
     read_pipeline(input_ch)
    }
}


// processed_reads_ch.subscribe{ println it }
