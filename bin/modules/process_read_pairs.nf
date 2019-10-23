//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process process_read_pairs {
  tag "$sampleID"
  beforeScript "module load gcc conda2"
  publishDir "$params.out_dir/split_reads", mode: "copy"

  input:
  tuple sampleID, file(reads)

  output:
  tuple sampleID, file('*_paired.fastq'), file('*_unpaired.fastq')

  script:
  """
  python $workflow.projectDir/bin/python/process_read_pairs.py \
  -f ${reads} \
  -o ${sampleID} \
  -l ${sampleID}_process_read_pairs.log
  """
}
