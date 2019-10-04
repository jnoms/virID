//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process bwa_mem_contigs {
  tag "$sampleID"
  beforeScript "module load gcc bwa samtools"

  publishDir "$params.out_dir/read_mapping", mode: "copy"

  input:
  set sampleID, contigs, paired_reads, unpaired_reads

  output:
  set sampleID, file("*_mapped.counts"), file("*_mapped.bam*"), file("*_unmapped.bam*")

  script:
  """
  $workflow.projectDir/bin/bash/map_contigs.sh \
  -p $paired_reads \
  -u $unpaired_reads \
  -c $contigs \
  -i ${sampleID}_index \
  -m ${sampleID}_mapped.bam \
  -U ${sampleID}_unmapped.bam \
  -o ${sampleID}_mapped.counts \
  -t ${task.cpus} \
  -l ${sampleID}_bwa.log \
  -e temp
  """
}
