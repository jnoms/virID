//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.source = "diamond"

//============================================================================//
// Define process
//============================================================================//
process get_counts {
  tag "$sampleID"
  publishDir "${params.out_dir}/${params.source}", mode: "copy"
  beforeScript "module load gcc conda2"

  input:
  tuple sampleID, assignment_file

  output:
  tuple sampleID, file("*tsv")

  script:
  """
  python $workflow.projectDir/bin/python/get_counts.py \
  -i ${assignment_file} \
  -o ${sampleID}_{$params.source}_counts.tsv \
  -l ${sampleID}_{$params.source}_counts.log
  """
}
