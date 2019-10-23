//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process fastq_to_fasta {
  tag "$sampleID"
  publishDir "$params.out_dir/fastq_to_fasta", mode: "copy"
  beforeScript "module load gcc conda2"

  input:
  tuple sampleID, fastq

  output:
  tuple sampleID, file("*fasta")

  script:
  """
  seqtk seq -a ${fastq} > ${sampleID}.fasta
  """
}
