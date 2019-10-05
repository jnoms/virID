//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"


//============================================================================//
// Define process
//============================================================================//
process bwa_mem_contigs {
  tag "$sampleID"
  beforeScript "module load gcc conda2"
  publishDir "$params.out_dir/read_mapping", mode: "copy"

  input:
  set sampleID, contigs, paired_reads, unpaired_reads

  output:
  set sampleID, file("*_mapped.counts"), file("*_cov"), file("*.bam")

  script:
  """
  $workflow.projectDir/bin/bash/map_contigs.sh \
  -p $paired_reads \
  -u $unpaired_reads \
  -c $contigs \
  -o ${sampleID}_mapped.counts \
  -v ${sampleID}_cov \
  -b ${sampleID}.bam \
  -t ${task.cpus} \
  -l ${sampleID}_bwa.log \
  -e temp
  """
}
