#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2


module load gcc
module load samtools
module load bwa

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#

usage() {
        echo "
        This purpose of this file is to map input reads to assembled contigs.
        Input is a fasta containing contigs and fastq(s) containing reads.
        Resultant counts file is tab-delimited with a structure of
        <contig_name> <counts>.

        USAGE: $0

        Required:
        -p <PAIRED_READS>        Path to the paired reads. PAIRED_READS,
                                 UNPAIRED_READS, or both are required.
        -u <UNPAIRED_READS>      Path to unpaired reads. PAIRED_READS,
                                 UNPAIRED_READS, or both are required.
        -c <CONTIGS>             Path to the contigs.
        -o <OUTPUT_COUNTS>       Path to the outfile containing counts.

        Optional:
        -i <INDEX_PATH>          Path to the contig index.
                                  Default: ./index
        -m <MAPPED_OUTPUT_BAM>   Path to the output BAM containing mapped reads
                                  Default: ./mapped.bam
        -u <UNMAPPED_OUTPUT_BAM> Path to the output BAM containing unmapped reads
                                  Default: ./unmapped.bam
        -t <THREADS>             Number of threads
                                  Default: 1
        -l <LOG_FILE>            Path to the log file.
                                  Default: No log file
        -e <TEMP_DIR>            Path to the temprorary directory
                                  Default: ./mapping

        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts p:u:c:i:m:U:o:t:l:e: option ; do
        case "${option}"
        in
                p) PAIRED_READS=${OPTARG};;
                u) UNPAIRED_READS=${OPTARG};;
                c) CONTIGS=${OPTARG};;
                i) INDEX_PATH=${OPTARG};;
                m) MAPPED_OUTPUT_BAM=${OPTARG};;
                U) UNMAPPED_OUTPUT_BAM=${OPTARG};;
                o) OUTPUT_COUNTS=${OPTARG};;
                t) THREADS=${OPTARG};;
                l) LOG_FILE=${OPTARG};;
                e) TEMP_DIR=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Defaults
#------------------------------------------------------------------------------#
INDEX_PATH=${INDEX_PATH:-./index}
MAPPED_OUTPUT_BAM=${MAPPED_OUTPUT_BAM:-./mapped.bam}
UNMAPPED_OUTPUT_BAM=${UNMAPPED_OUTPUT_BAM:-./unmapped.bam}
THREADS=${THREADS:-1}
TEMP_DIR=${TEMP_DIR:-./mapping}

#------------------------------------------------------------------------------#
#Defining functions
#------------------------------------------------------------------------------#
write_log() {
  #The purpose of this function is to echo a message to stdout, and write that message to a log file with the date. This function works as follows:
  #write_log $1="message" $2=<log_file> $3="error"(optional)
  local MESSAGE=$1
  local LOG_FILE=$2
  local ERROR_STATUS=$3

  # Print message to screen
  echo "$MESSAGE"
  date

  # If there is no log file, we're done here
  if [[ $LOG_FILE == '' ]] ; then
    return
  fi

  # Print message to log file
  echo "$MESSAGE" >> $LOG_FILE
  date >> $LOG_FILE

  # If this is an error, write to error log file.
  if [[ $ERROR_STATUS == "error" ]] ; then
    echo "$MESSAGE" >> $LOG_FILE'_error'
    date >> $LOG_FILE'_error'
  fi
}

check_if_file_exists() {
  #This function checks if input files exist. If they don't, it prints an error message and exits script.

  #Input variables (must exist):
  # $1 - A space-delimited list of files to check for existance, ex "path_to_file1 path_to_file2 path_to_file3"
  # $LOG_FILE - path to the log file. Assumed that the error log is $LOG_FILE'_error'

  #OUTPUT: Will write message to $LOG_FILE and $LOG_FILE'_error' detailing which file doesn't exist, and will exit script.

  FILE_LIST=$1
  for FILE in $FILE_LIST ; do
    if [ ! -s $FILE ] ; then
      echo "Cannot find the file $FILE. Exiting script."
      echo "Cannot find the file $FILE. Exiting script." >> $LOG_FILE
      echo "Cannot find the file $FILE. Exiting script." >> $LOG_FILE'_error'
      exit 1
    fi
  done
}

paired_or_unpaired() {
  #This file checks if the paired and/or unpaired files exist. If neither exists, exists script.

  #Input variables (must exist):
  # $PAIRED_READS - path to the file containing paired reads
  # $UNPAIRED_READS - path to the file containing unpaired reads
  # $LOG_FILE - path to the log file.
  # $SAMPLE - should be a string dictating the sample name.

  # Function depedencies:
  # write_log

  #OUTPUT: Setting the variable $PAIRED_OR_UNPAIRED as either "PAIRED", "UNPAIRED", or "BOTH"

  if [[ -s $PAIRED_READS ]] && [[ -s $UNPAIRED_READS ]] ; then
    export PAIRED_OR_UNPAIRED="BOTH"
  elif [[ -s $PAIRED_READS ]] && [[ ! -s $UNPAIRED_READS ]] ; then
    export PAIRED_OR_UNPAIRED="PAIRED"
  elif [[ ! -s $PAIRED_READS ]] && [[ -s $UNPAIRED_READS ]] ; then
    export PAIRED_OR_UNPAIRED="UNPAIRED"
  elif [[ ! -s $PAIRED_READS ]] && [[ ! -s $UNPAIRED_READS ]] ; then
    echo "Can't find either a paired or unpaired file!. Exiting."
    echo "Can't find either a paired or unpaired file!. Exiting." >> $LOG_FILE
    echo "Can't find either a paired or unpaired file!. Exiting." >> $LOG_FILE'_error'
    exit 1
  fi

  write_log "map_contigs.bash, $SAMPLE: Pairing set to $PAIRED_OR_UNPAIRED" $LOG_FILE
}

make_index() {
  #Make BWA-mem index from input contigs. Required variables:

  #Input variables (must exist):
  # $CONTIGS - path to the file containing contigs.
  # $INDEX_PATH - desired path to the index. Format should be: path/you/want/index_name
  # $SAMPLE - should be a string dictating the sample name.
  # $LOG_FILE - path to the log file.

  #Function dependencies:
  # write_log
  # check_if_file_exists

  #OUTPUT: A BWA index at $INDEX_PATH

  write_log "map_contigs.bash, $SAMPLE: Making BWA-MEM index out of $CONTIGS." $LOG_FILE

  local SAMPLE=$( basename $CONTIGS )
  local SAMPLE=${SAMPLE%.fasta}

  local INDEX_DIR=$(dirname $INDEX_PATH)
  local INDEX_NAME=$(basename $INDEX_PATH)

  #Make the index directory
  mkdir -p $INDEX_DIR

  #Make the index
  bwa index -p $INDEX_NAME $CONTIGS

  #move the index files to the index directory - working index is $INDEX_PATH
  mv $INDEX_NAME* $INDEX_DIR

  #Make sure the sample index exists
  check_if_file_exists $INDEX_PATH'.amb'

  write_log "map_contigs.bash, $SAMPLE: Finished making BWA-MEM index out of $CONTIGS." $LOG_FILE
}

run_BWA_mem_paired() {
  #Runs BWA-MEM on paired, interleaved input reads.

  #Input variables (must exist):
  # $INDEX_PATH - path to the index. Format should be: path/to/the/index_name
  # $PAIRED_READS - path to the paired reads
  # $LOG_FILE - path to the log file.
  # $1 - path to output bam

  #Function dependencies:
  # check_if_file_exists

  check_if_file_exists $PAIRED_READS

  bwa mem \
  -t $THREADS \
  -p \
  -k 15 \
  $INDEX_PATH \
  $PAIRED_READS | samtools view -bh - > $1
}

run_BWA_mem_unpaired() {
  #Runs BWA-MEM on unpaired input reads.

  #Input variables (must exist):
  # $INDEX_PATH - path to the index. Format should be: path/to/the/index_name
  # $UNPAIRED_READS - path to the paired reads
  # $1 - path to output bam

  #Function dependencies:
  # check_if_file_exists

  check_if_file_exists $UNPAIRED_READS

  bwa mem \
  -t $THREADS \
  -k 15 \
  $INDEX_PATH \
  $UNPAIRED_READS | samtools view -bh - > $1
}

run_BWA_mem() {
  #Runs BWA mem on input contigs and reads.

  #Input: The following variables:
  # $PAIRED_OR_UNPAIRED - this variable should be set to "PAIRED", "UNPAIRED", or "BOTH"
  # $SAMPLE - name of the sample
  # $TEMP_DIR - path to the temporary directory
  # $1 - path to the output BAM (which includes mapped and unmapped reads)

  #Function dependencies:
  # write_log
  # check_if_file_exists

  #Output: Bam, including mapped and unmapped reads, will be output to the path specified in $1

  write_log "map_contigs.bash, $SAMPLE: Running BWA-mem to map reads to contigs." $LOG_FILE

  if [[ $PAIRED_OR_UNPAIRED == "BOTH" ]] ; then
    run_BWA_mem_paired $TEMP_DIR/$SAMPLE'_paired_temp.bam'
    run_BWA_mem_unpaired $TEMP_DIR/$SAMPLE'_unpaired_temp.bam'

    samtools merge $1 $TEMP_DIR/$SAMPLE'_paired_temp.bam' $TEMP_DIR/$SAMPLE'_unpaired_temp.bam'

  elif [[ $PAIRED_OR_UNPAIRED == "PAIRED" ]] ; then
    run_BWA_mem_paired $1
  elif [[ $PAIRED_OR_UNPAIRED == "UNPAIRED" ]] ; then
    run_BWA_mem_unpaired $1
  fi

  check_if_file_exists $1

  write_log "map_contigs.bash, $SAMPLE: Finished running BWA-mem to map reads to contigs." $LOG_FILE
}

split_mapped-unmapped_and_sort_and_count() {
  #Isolates the mapped reads from a sam using samtools, and reports back the # of reads that are mapped and unmapped.
  #Also does samtools sort on each. It outputs a count for the mapped contigs.

  #Input: The following variables:
  # $LOG_FILE - Path to the log file
  # $TEMP_DIR - path to the temporary directory
  # $SAMPLE - name of the sample
  # $1 - Path to the input bam
  # $2 - Path to the output bam containing only mapped reads
  # $3 - Path to the output bam containing only unmapped reads
  # $4 - Path to the output counts file

  #Output: Sam only containing mapped reads at the path determined by $2

  local INPUT_BAM=$1
  local MAPPED_OUTPUT_BAM=$2
  local UNMAPPED_OUTPUT_BAM=$3
  local OUTPUT_COUNTS=$4

  #Determine the number of mapped and unmapped reads
  echo "Counting number of unmapped reads in $INPUT_BAM" >> $LOG_FILE
  samtools view -c -f 4 $INPUT_BAM >> $LOG_FILE

  echo "Counting number of mapped reads in $INPUT_BAM" >> $LOG_FILE
  samtools view -c -F 4 $INPUT_BAM >> $LOG_FILE

  #Taking only unmapped reads
  samtools view -h -b -f 4 $INPUT_BAM > $TEMP_DIR/$SAMPLE'_unmapped_unsorted'
  samtools sort $TEMP_DIR/$SAMPLE'_unmapped_unsorted' > $UNMAPPED_OUTPUT_BAM
  samtools index $UNMAPPED_OUTPUT_BAM

  #Taking only mapped reads
  samtools view -h -b -F 4 $INPUT_BAM > $TEMP_DIR/$SAMPLE'_mapped_unsorted'
  samtools sort $TEMP_DIR/$SAMPLE'_mapped_unsorted' > $MAPPED_OUTPUT_BAM
  samtools index $MAPPED_OUTPUT_BAM

  #Getting mapped read counts
  samtools idxstats $MAPPED_OUTPUT_BAM | cut -f1,3 > $OUTPUT_COUNTS

}

get_contig_counts() {
  #Gets number of reads that map to each contig. Input is a .sam, and output is a file with
  #the following structure: <contig_name> <count>

  #Input: The following variables
  # $LOG_FILE - Path to the log file
  # $SAMPLE - The sample name to be used for output naming
  # $1 - Path to the input bam
  # $2 - Path to the output counts file

  #Function dependencies:
  # write_log
  # check_if_file_exists

  #Output: A file at the path specified by $2 containing contig names and read counts.

  write_log "map_contigs.bash, $SAMPLE: Determining number of mapped reads for each contig." $LOG_FILE

  local INPUT_BAM=$1
  local OUTPUT_COUNTS=$2

  local BAM_DIR=$(dirname $INPUT_BAM)

  samtools sort $INPUT_BAM > $BAM_DIR/$SAMPLE'_sorted.bam'
  samtools index $BAM_DIR/$SAMPLE'_sorted.bam'
  samtools idxstats $BAM_DIR/$SAMPLE'_sorted.bam' | cut -f1,3 > $OUTPUT_COUNTS
  check_if_file_exists $OUTPUT_COUNTS

}


#------------------------------------------------------------------------------#
#Executing script
#------------------------------------------------------------------------------#
#Find sample name
SAMPLE=$(basename ${CONTIGS%.fasta})

#Make required directories
mkdir -p $(dirname $MAPPED_OUTPUT_BAM)
mkdir -p $(dirname $UNMAPPED_OUTPUT_BAM)
mkdir -p $(dirname $LOG_FILE)
mkdir -p $(dirname $OUTPUT_COUNTS)
mkdir -p $TEMP_DIR

#Updating log
write_log "map_contigs.bash, $SAMPLE: Starting map_contigs.bash." $LOG_FILE

#Check that input files exist
check_if_file_exists "$CONTIGS"

#Main
paired_or_unpaired
make_index
run_BWA_mem $TEMP_DIR/$SAMPLE'_mapped_and_unmapped_reads.bam'
split_mapped-unmapped_and_sort_and_count $TEMP_DIR/$SAMPLE'_mapped_and_unmapped_reads.bam' $MAPPED_OUTPUT_BAM $UNMAPPED_OUTPUT_BAM $OUTPUT_COUNTS
#get_contig_counts $MAPPED_OUTPUT_BAM $OUTPUT_COUNTS

#Make sure the output files exists
check_if_file_exists "$MAPPED_OUTPUT_BAM"

#updating log
write_log "map_contigs.bash, $SAMPLE: Completed map_contigs.bash." $LOG_FILE
