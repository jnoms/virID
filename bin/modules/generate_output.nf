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
  set sampleID, blast_file, diamond_file, contigs, mapped_counts, mapped_coverage, mapped_bam, contaminant

  output:
  set sampleID, file("*tsv")

  script:
  """
  cut -f1 ${contaminant} > contaminant_list.txt

  # If contaminant file has contents, subtract those sequences before generating
  # counts.
  if [[ -s contaminant_list.txt ]] ; then
    grep -v -f contaminant_list.txt ${blast_file} > blast_contaminants_removed.txt
    grep -v -f contaminant_list.txt ${diamond_file} > diamond_contaminants_removed.txt
  else
    cp ${blast_file} blast_contaminants_removed.txt
    cp ${diamond_file} diamond_contaminants_removed.txt
  fi

  # Generate counts files - blast
  python $workflow.projectDir/bin/python/get_counts.py \
  -i blast_contaminants_removed.txt \
  -o ${sampleID}_blast_counts.tsv \
  -l ${sampleID}_get_counts.log \
  -c ${mapped_counts}

  # Generate counts files - diamond
  python $workflow.projectDir/bin/python/get_counts.py \
  -i diamond_contaminants_removed.txt \
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
  -v ${mapped_coverage} \
  -m contaminant_list.txt \
  -l ${sampleID}_merge.log
  """
}
