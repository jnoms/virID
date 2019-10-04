//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process generate_output {
  tag "$sampleID"
  publishDir "$params.out_dir/results", mode: "copy"
  beforeScript "module load gcc conda2 ; source activate conda_py36"

  input:
  set sampleID, blast_file, diamond_file, contigs, mapped_counts, mapped_bam, unmapped_bam

  output:
  set sampleID, file("*tsv")

  script:
  """
  # Generate counts files - blast
  python $workflow.projectDir/bin/python/get_counts.py \
  -i ${blast_file} \
  -o ${sampleID}_blast_counts.tsv \
  -l ${sampleID}_get_counts.log \
  -c ${mapped_counts}

  # Generate counts files - diamond
  python $workflow.projectDir/bin/python/get_counts.py \
  -i ${diamond_file} \
  -o ${sampleID}_diamond_counts.tsv \
  -l ${sampleID}_get_counts.log \
  -c ${mapped_counts}

  # Generate merged output file
  python $workflow.projectDir/bin/python/merge_and_postprocess.py \
  -t ${blast_file} \
  -T ${diamond_file} \
  -f ${contigs} \
  -o ${sampleID}_merged.tsv \
  -c ${mapped_counts} \
  -l ${sampleID}_merge.log
  """
}
