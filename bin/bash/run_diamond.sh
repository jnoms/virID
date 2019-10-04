#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2

module load gcc

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#

usage() {
        echo "
        This purpose of this file is to use DIAMOND to map a file containing
        query sequences against a pre-determined protein database. DIAMOND will
        automatically detect and use # of threads available.

        USAGE: $0

        Required:
        -d <DATABASE>    Path to the diamond database.
        -q <QUERY>       Path to the query file containing nucleotide sequences.
        -o <OUTPUT_PATH> Path to the output file.

        Optional:
        -m <MEMORY>      Memory available in GB. DIAMOND blocksize is set to
                          Default: 10
        -t <TEMP_DIR>    Path to the temporary directory
                          Default: ./diamond_temp
        -e <EVALUE>      Max evalue.
                          Default: 10
        -f <OUT_FORMAT>  Output format.
                          Default: '6 qseqid stitle sseqid staxids evalue
                          bitscore pident length'
        -l <LOG_FILE>   Path to the log file
                           Default: No log file.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts d:q:o:m:t:e:f:l: option ; do
        case "${option}"
        in
                d) DATABASE=${OPTARG};;
                q) QUERY=${OPTARG};;
                o) OUTPUT_PATH=${OPTARG};;
                m) MEMORY=${OPTARG};;
                t) TEMP_DIR=${OPTARG};;
                e) EVALUE=${OPTARG};;
                f) OUT_FORMAT=${OPTARG};;
                l) LOG_FILE=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Setting defaults
#------------------------------------------------------------------------------#
MEMORY=${MEMORY:-10}
TEMP_DIR=${TEMP_DIR:-./diamond_temp}
EVALUE=${EVALUE:-10}
OUT_FORMAT=${OUT_FORMAT:-"6 qseqid stitle sseqid staxids evalue bitscore pident length"}

# Calculate DIAMOND block size
if [[ $MEMORY < 10 ]] ; then
  BLOCK=1
else
  BLOCK=$((MEMORY/10))
fi

#------------------------------------------------------------------------------#
#Defining functions
#------------------------------------------------------------------------------#
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

run_diamond() {
  #Make sure all necesssary directories have been made
  mkdir -p $( dirname $OUTPUT_PATH )
  mkdir -p $TEMP_DIR

  #Get sample name
  SAMPLE=$( basename $QUERY )

  #Write to log
  write_log "run_diamond.bash, $SAMPLE: Starting run_diamond.bash script. " $LOG_FILE

  #Make sure inputs exist exists
  if [ ! -f $DATABASE ] ; then
    write_log "run_diamond.bash, $SAMPLE: Cannot find the input database, $DATABASE. Exiting." $LOG_FILE "error"
    exit 1
  elif [ ! -f $QUERY ] ; then
    write_log "run_diamond.bash, $SAMPLE: Cannot find the query file, $QUERY. Exiting." $LOG_FILE "error"
    exit 1
  fi

  #Running DIAMOND
  diamond blastx \
  -d $DATABASE \
  -q $QUERY \
  --more-sensitive \
  -o $OUTPUT_PATH \
  --tmpdir $TEMP_DIR \
  --evalue $EVALUE \
  --outfmt $OUT_FORMAT \
  --index-chunks 1 \
  --top 1 \
  --block-size $BLOCK

  write_log "run_diamond.bash, $SAMPLE: Finished run_diamond.bash script." $LOG_FILE

  #Check that the output file exists. It no matches it should exist, though it will be empty.
  if [ ! -f $OUTPUT_PATH ] ; then
    write_log "run_diamond.bash, $SAMPLE: CANNOT find DIAMOND output file. Exiting script." $LOG_FILE "error"
    exit 1
  fi
}

#------------------------------------------------------------------------------#
#Calling functions
#------------------------------------------------------------------------------#
run_diamond
