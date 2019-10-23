//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.log_file = "${workflow.launchDir}/${params.out_dir}/reports/virID.log"

//============================================================================//
// Define process
//============================================================================//
process bwa_mem_contigs {
  tag "$sampleID"
  beforeScript "module load gcc conda2"
  publishDir "$params.out_dir/read_mapping", mode: "copy"

  input:
  tuple sampleID, contigs, paired_reads, unpaired_reads

  output:
  tuple sampleID, file("*_mapped.counts"), file("*_cov"), file("*.bam")

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
  -l ${params.log_file} \
  -e temp \
  -s ${sampleID}
  """
}
