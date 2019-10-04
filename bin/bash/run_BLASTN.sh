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
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts d:q:o:t:e:f:l:b:m:s:r:i: option ; do
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

###############################
#Defining functions
###############################
write_log() {
  # The purpose of this function is to echo a message to stdout, and write that
  # message to a log file with the date. This function works as follows:
  # write_log $1="message" $2=<log_file> $3="error"(optional)
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
  mkdir -p $(dirname $LOG_FILE)
  echo "$MESSAGE" >> $LOG_FILE
  date >> $LOG_FILE

  # If this is an error, write to error log file.
  if [[ $ERROR_STATUS == "error" ]] ; then
    echo "$MESSAGE" >> $LOG_FILE'_error'
    date >> $LOG_FILE'_error'
  fi
}

run_BLASTN() {
  #Get sample name
  SAMPLE=$( basename $QUERY )

  #Make sure that the necessary directories have been made
  mkdir -p $( dirname $OUTPUT_FILE )

  #Write to log
  write_log "run_BLASTN.bash, $SAMPLE: Starting run_BLASTN script." $LOG_FILE

  #Check that inputs exist
  if [[ ! -s $QUERY ]] ; then
    write_log "run_BLASTN.bash, $SAMPLE: Cannot find QUERY file, which is set to $QUERY. Exiting." $LOG_FILE "error"
    exit 1
  fi

  #Call BLASTN
  if [[ $RESTRICT_TO_TAXIDS == '' ]] && [[ $IGNORE_TAXIDS == '' ]] ; then
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
  elif [[ $RESTRICT_TO_TAXIDS != '' ]] && [[ $IGNORE_TAXIDS == '' ]] ; then
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
  elif [[ $RESTRICT_TO_TAXIDS == '' ]] && [[ $IGNORE_TAXIDS != '' ]] ; then
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
  fi

  #Check that the output file exists. It no matches it should exist, though it will be empty.
  if [[ ! -f $OUTPUT_FILE ]] ; then
    write_log "run_BLASTN.bash, $SAMPLE: Cannot find output file, suggesting something went wrong. Exiting." $LOG_FILE "error"
    exit 1
  fi

  #Write to log
  write_log "run_BLASTN.bash, $SAMPLE: Finished run_BLASTN script." $LOG_FILE
}

###############################
#Calling functions
###############################
run_BLASTN
