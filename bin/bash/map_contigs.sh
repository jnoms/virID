#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2

module load gcc
module load samtools
module load bwa

set -e

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        This purpose of this file is to map input reads to assembled contigs.
        Input is a fasta containing contigs and fastq(s) containing reads.

        Output files are:
        1) a sorted bam containing the alignments,
        2) a file detailing read counts - <ID> <count> (tsv, no header)
        3) a file detailing coverage info - <#ID> <Avg_fold> <Covered_percent>
            (space-delimited, with header)

        USAGE: $0

        Required:
        -p <PAIRED_READS>        Path to the paired reads. PAIRED_READS,
                                 UNPAIRED_READS, or both are required.
        -u <UNPAIRED_READS>      Path to unpaired reads. PAIRED_READS,
                                 UNPAIRED_READS, or both are required.
        -c <CONTIGS>             Path to the contigs.
        -o <OUTPUT_COUNTS>       Path to the outfile containing counts.
                                 Format: ID counts <no header>
        -v <OUTPUT_COV>          Path to output coverage file (space-delimited).
                                 Format: #ID Avg_fold Covered_percent
        -b <OUTPUT_BAM>          Path to the output bam.

        Optional:
        -t <THREADS>             Number of threads
                                  Default: 1
        -l <LOG_FILE>            Path to the log file.
                                  Default: No log file
        -e <TEMP_DIR>            Path to the temprorary directory
                                  Default: ./mapping
        -s <SAMPLE_ID>           Name of the sample being processed.
                                  Default: Basename of PAIRED_READS

        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts p:u:c:o:v:b:t:l:e:s: option ; do
        case "${option}"
        in
                p) PAIRED_READS=${OPTARG};;
                u) UNPAIRED_READS=${OPTARG};;
                c) CONTIGS=${OPTARG};;
                o) OUTPUT_COUNTS=${OPTARG};;
                v) OUTPUT_COV=${OPTARG};;
                b) OUTPUT_BAM=${OPTARG};;

                t) THREADS=${OPTARG};;
                l) LOG_FILE=${OPTARG};;
                e) TEMP_DIR=${OPTARG};;
                s) SAMPLE_ID=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Defaults
#------------------------------------------------------------------------------#
THREADS=${THREADS:-1}
TEMP_DIR=${TEMP_DIR:-./mapping}
LOG_FILE=${LOG_FILE:-$TEMP_DIR/log.txt}
SAMPLE_ID=${SAMPLE_ID:-$(basename $PAIRED_READS)}

#------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------#
write_log() {
  # The purpose of this function is to append to a log, which is CSV output.
  # Structure is sampleID, file_name/process, message, error(0 or 1). Log file
  # is optional - if not input, will print to screen.

  # Usage:
  # write_log <sampleID> <message> <error - 0 or 1> <log_file>
  sampleID=$1
  FILE_NAME=$(basename $0)
  message=$2
  error=$3
  log_file=$4

  OUTPUT="${sampleID},${FILE_NAME},${message},${error}"
  echo $OUTPUT
  date

  # If no log file, leave
  if [[ $log_file == '' ]] ; then
    return
  fi

  # Write to log
  mkdir -p $(dirname $log_file)
  echo $OUTPUT >> $log_file
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

      write_log \
      $SAMPLE_ID \
      "Cannot find the file $FILE. Exiting script." \
      1 \
      $LOG_FILE

      exit 1
    fi
  done
}

run_bwa() {
  READS=$1
  OUT_BAM=$2

  # If reads are empty, put out an empty file
  if [[ ! -s $READS ]] ; then
    touch $OUT_BAM
    return
  fi

  bwa mem \
  -t $THREADS \
  -p \
  -k 15 \
  $TEMP_DIR/index \
  $READS | samtools view -bh - | samtools sort > $OUT_BAM
}

#------------------------------------------------------------------------------#
# Execute script
#------------------------------------------------------------------------------#

write_log \
$SAMPLE_ID \
"Starting script." \
0

# Set up dirs
#------------------------------------------------------------------------------#
mkdir -p $TEMP_DIR
mkdir -p $(dirname $LOG_FILE)
mkdir -p $(dirname $OUTPUT_BAM)
mkdir -p $(dirname $OUTPUT_COV)
mkdir -p $(dirname $OUTPUT_COUNTS)

# Make index
#------------------------------------------------------------------------------#
bwa index -p $TEMP_DIR/index $CONTIGS
check_if_file_exists $TEMP_DIR/index.amb

# Run BWA
#------------------------------------------------------------------------------#
run_bwa $PAIRED_READS $TEMP_DIR/paired.bam
run_bwa $UNPAIRED_READS $TEMP_DIR/unpaired.bam

# If there are both files, merge them; else use one
if [[ -s $TEMP_DIR/paired.bam ]] &&  [[ -s $TEMP_DIR/unpaired.bam ]] ; then
  samtools merge $TEMP_DIR/mapped.bam $TEMP_DIR/paired.bam $TEMP_DIR/unpaired.bam
elif [[ -s $TEMP_DIR/paired.bam ]] &&  [[ ! -s $TEMP_DIR/unpaired.bam ]] ; then
  cp $TEMP_DIR/paired.bam $TEMP_DIR/mapped.bam
elif [[ ! -s $TEMP_DIR/paired.bam ]] &&  [[ -s $TEMP_DIR/unpaired.bam ]] ; then
  cp $TEMP_DIR/unpaired.bam $TEMP_DIR/mapped.bam
fi
check_if_file_exists $TEMP_DIR/mapped.bam

# Sort, count, and get cov
#------------------------------------------------------------------------------#
samtools sort $TEMP_DIR/mapped.bam > $OUTPUT_BAM
samtools index $OUTPUT_BAM

samtools idxstats $OUTPUT_BAM | cut -f1,3 > $OUTPUT_COUNTS
pileup.sh in=${OUTPUT_BAM} out=$TEMP_DIR/coverage.txt
awk '{print $1,$2,$5}' $TEMP_DIR/coverage.txt > $OUTPUT_COV

write_log \
$SAMPLE_ID \
"Finished." \
0
