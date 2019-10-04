//============================================================================//
// Default params
//============================================================================//
params.within_percent_of_top_score = 1
params.out_dir = "output"
params.column_names = "query_ID seq_title seq_ID taxonID evalue bitscore pident length"
params.source = "diamond"
params.taxid_blacklist = "$workflow.projectDir/resources/2019-08-09_blacklist.tsv"
params.taxonomy_column = "taxonID"
params.score_column = "bitscore"

//============================================================================//
// Define process
//============================================================================//
process blast_to_LCA {
  tag "$sampleID"
  publishDir "${params.out_dir}/${params.source}", mode: "copy"
  beforeScript "module load gcc conda2 ; source activate conda_py36"

  input:
  set sampleID, assignment_file

  output:
  set sampleID, file("*.tsv")

  script:
  """
  python $workflow.projectDir/bin/python/get_LCA.py \
  -i ${assignment_file} \
  -o ${sampleID}_${params.source}.tsv \
  -c "${params.column_names}" \
  -t ${params.taxonomy_column} \
  -s ${params.score_column} \
  -b ${params.taxid_blacklist} \
  -w ${params.within_percent_of_top_score} \
  -l ${sampleID}_conversion.log
  """
}
