#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2

set -e

#------------------------------------------------------------------------------#
#Defining usage and setting input
#------------------------------------------------------------------------------#

usage() {
        echo "
        This purpose of this file is to run blast on an input query file.

        USAGE: $0

        Required:
        -d <DATABASE>    Path to BLAST database.
        -q <QUERY>       Path to the query file
        -o <OUTPUT_FILE> Path to the BLASTN output file

        Optional:
        -t <THREADS>     Number of threads available              Default: 1
        -e <EVALUE>      Maximum evalue.                          Default: 10
        -f <OUT_FORMAT>  BLAST output format.                     Default:
                              '6 qseqid stitle sseqid staxid evalue bitscore pident length'
        -l <LOG_FILE>    Path to and name of the log file         Default: No log file.
        -b <BLAST_TYPE>  BLAST task. Options are megablast, dc-megablast, and
                         blastn.                                  Default: megablast
        -m <MAX_HSPS>    Maximum number of alignments per hit.    Default: 1
        -s <MAX_TARGETS> Maximum number of hits per query.        Default: 30
        -r <RESTRICT_TO_TAXIDS> Only search against this taxonID. Default: None
        -i <IGNORE_TAXIDS>      Ignore these taxonIDs.            Default: None
        -n <SAMPLE_ID>   Name of the sample being processed.      Default: Basename of QUERY
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts d:q:o:t:e:f:l:b:m:s:r:i:n: option ; do
        case "${option}"
        in
                d) DATABASE=${OPTARG};;
                q) QUERY=${OPTARG};;
                o) OUTPUT_FILE=${OPTARG};;
                t) THREADS=${OPTARG};;
                e) EVALUE=${OPTARG};;
                f) OUT_FORMAT=${OPTARG};;
                l) LOG_FILE=${OPTARG};;
                b) BLAST_TYPE=${OPTARG};;
                m) MAX_HSPS=${OPTARG};;
                s) MAX_TARGETS=${OPTARG};;
                r) RESTRICT_TO_TAXIDS=${OPTARG};;
                i) IGNORE_TAXIDS=${OPTARG};;
                n) SAMPLE_ID=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Setting defaults
#------------------------------------------------------------------------------#
THREADS=${THREADS:-1}
EVALUE=${EVALUE:-10}
OUT_FORMAT=${OUT_FORMAT:-'6 qseqid stitle sseqid staxid evalue bitscore pident length'}
BLAST_TYPE=${BLAST_TYPE:-megablast}
MAX_HSPS=${MAX_HSPS:-1}
MAX_TARGETS=${MAX_TARGETS:-30}
RESTRICT_TO_TAXIDS=${RESTRICT_TO_TAXIDS:-no}
IGNORE_TAXIDS=${IGNORE_TAXIDS:-no}
SAMPLE_ID=${SAMPLE_ID:-$(basename $QUERY)}


###############################
#Defining functions
###############################
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

run_BLASTN() {
  #Get sample name
  SAMPLE=$( basename $QUERY )

  #Make sure that the necessary directories have been made
  mkdir -p $( dirname $OUTPUT_FILE )

  #Write to log
  write_log \
  $SAMPLE_ID \
  "Starting script at $(date)." \
  0

  #Check that inputs exist
  if [[ ! -s $QUERY ]] ; then
    write_log \
    $SAMPLE_ID \
    "Cannot find QUERY file, which is set to $QUERY. Exiting." \
    1 \
    $LOG_FILE

    exit 1
  fi

  #Call BLASTN
  if [[ $RESTRICT_TO_TAXIDS == 'no' ]] && [[ $IGNORE_TAXIDS == 'no' ]] ; then
    blastn \
    -query $QUERY \
    -db $DATABASE \
    -out $OUTPUT_FILE \
    -outfmt "$OUT_FORMAT" \
    -num_threads $THREADS \
    -evalue $EVALUE \
    -task $BLAST_TYPE \
    -max_hsps $MAX_HSPS \
    -max_target_seqs $MAX_TARGETS
  elif [[ $RESTRICT_TO_TAXIDS != 'no' ]] && [[ $IGNORE_TAXIDS == 'no' ]] ; then
    blastn \
    -query $QUERY \
    -db $DATABASE \
    -out $OUTPUT_FILE \
    -outfmt "$OUT_FORMAT" \
    -num_threads $THREADS \
    -evalue $EVALUE \
    -task $BLAST_TYPE \
    -max_hsps $MAX_HSPS \
    -max_target_seqs $MAX_TARGETS \
    -taxids $RESTRICT_TO_TAXIDS
  elif [[ $RESTRICT_TO_TAXIDS == 'no' ]] && [[ $IGNORE_TAXIDS != 'no' ]] ; then
    blastn \
    -query $QUERY \
    -db $DATABASE \
    -out $OUTPUT_FILE \
    -outfmt "$OUT_FORMAT" \
    -num_threads $THREADS \
    -evalue $EVALUE \
    -task $BLAST_TYPE \
    -max_hsps $MAX_HSPS \
    -max_target_seqs $MAX_TARGETS \
    -negative_taxids $IGNORE_TAXIDS
  else
    write_log \
    $SAMPLE_ID \
    "Cannot specify both RESTRICT_TO_TAXIDS and IGNORE_TAXIDS." \
    1 \
    $LOG_FILE

    exit 1
  fi

  #Check that the output file exists. It no matches it should exist, though it will be empty.
  if [[ ! -f $OUTPUT_FILE ]] ; then
    write_log \
    $SAMPLE_ID \
    "Cannot find output file, suggesting something went wrong. Exiting." \
    1 \
    $LOG_FILE

    exit 1
  fi

  # If there is no contents in the OUTPUT_FILE don't quit, but report warning
  if [[ ! -s $OUTPUT_PATH ]] ; then
    write_log \
    $SAMPLE_ID \
    "WARN: Output file $OUTPUT_FILE doesn't have any contents." \
    1 \
    $LOG_FILE
  fi

  #Write to log
  write_log \
  $SAMPLE_ID \
  "Finished script at $(date)." \
  0
}

###############################
#Calling functions
###############################
run_BLASTN
