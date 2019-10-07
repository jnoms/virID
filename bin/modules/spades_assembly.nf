//============================================================================//
// Default params
//============================================================================//
params.spades_type = 'rna'
params.temp_dir = 'temp'
params.out_dir = 'output'
params.spades_min_length = 300
params.log_file = "${launchDir}/${params.out_dir}/results/virID.log"

//============================================================================//
// Define process
//============================================================================//
process spades_assembly {
  tag "$sampleID"
  publishDir "$params.out_dir/contigs", mode: "copy"
  beforeScript "module load gcc conda2"

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
  -l ${params.log_file} \
  -m ${task.memory.toGiga()} \
  -t ${task.cpus} \
  -e ${params.temp_dir} \
  -o ${sampleID}_contigs.fasta \
  -L ${params.spades_min_length} \
  -n ${sampleID}
  """
}
