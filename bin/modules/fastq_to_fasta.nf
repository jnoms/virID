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
  tuple sampleID, paired_fastq, unpaired_fastq

  output:
  tuple sampleID, file("*fasta")

  script:
  """
  # Concatenate the paired and unpaired files. The entire purpose of splitting
  # them is to make sure paired reads are appropriately marked with /1 and /2.
  cat ${paired_fastq} ${unpaired_fastq} > reads.fastq

  # Convert to fasta
  seqtk seq -a reads.fastq > ${sampleID}.fasta
  """
}
