//============================================================================//
// Default params
//============================================================================//
params.diamond_database = "ALL_SMALL80_NH"
params.diamond_outfmt = "6 qseqid stitle sseqid staxids evalue bitscore pident length"
params.temp_dir = "temp"
params.diamond_evalue = "10"
params.out_dir = "output"
params.log_file = "${workflow.launchDir}/${params.out_dir}/reports/virID.log"

//============================================================================//
// Define process
//============================================================================//
process diamond {
  tag "$sampleID"
  publishDir "$params.out_dir/diamond", mode: "copy"
  beforeScript "module load gcc conda2"

  input:
  tuple sampleID, sequences

  output:
  tuple sampleID, file("*_diamond.out")

  script:
  """
  $workflow.projectDir/bin/bash/run_diamond.sh \
  -d ${params.diamond_database} \
  -q ${sequences} \
  -o ${sampleID}_diamond.out \
  -m ${task.memory.toGiga()} \
  -t ${params.temp_dir} \
  -e ${params.diamond_evalue} \
  -f "${params.diamond_outfmt}" \
  -l ${params.log_file} \
  -s ${sampleID}
  """
}
