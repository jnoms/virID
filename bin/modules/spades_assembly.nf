//============================================================================//
// Default params
//============================================================================//
params.spades_type = 'rna'
params.temp_dir = 'temp'
params.out_dir = 'output'

//============================================================================//
// Define process
//============================================================================//
process spades_assembly {
  tag "$sampleID"
  publishDir "$params.out_dir/contigs", mode: "copy"

  input:
  set sampleID, paired_reads, unpaired_reads

  output:
  set sampleID, file('*contigs.fasta')

  script:
  """
  $workflow.projectDir/bin/bash/SPAdes.sh \
  -s ${params.spades_type} \
  -p ${paired_reads} \
  -u ${unpaired_reads} \
  -l ${sampleID}_SPAdes_assembly.log \
  -m ${task.memory.toGiga()} \
  -t ${task.cpus} \
  -e ${params.temp_dir} \
  -o ${sampleID}_contigs.fasta \
  -g yes
  """
}
