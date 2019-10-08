#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#

usage() {
        echo "
        This purpose of this file is to take an input .fastq and assemble it
        using SPAdes.

        USAGE: $0

        Required:
        -p  <PAIRED_FILE>      Path to file. Must specify PAIRED_FILE,
                               UNPAIRED_FILE, or both.
        -u  <UNPAIRED_FILE>    Path to file. Must specify PAIRED_FILE,
                               UNPAIRED_FILE, or both.
        -o  <OUTPUT_FILE_NAME> Path to and name of the output file. Recommend
                               something like SAMPLE'_contigs.fasta.'

        Optional:
        -s  <SPADES_TYPE>     meta, rna, or regular.   Default: rna
        -l  <LOG_FILE>        Log file.                Default: no log file
        -m  <MEMORY>          Number of GB of memory.  Default: 10
        -t  <THREADS>         Number of threads        Default: 1
        -e  <TEMP_DIR>        TEMP directory.          Default: ./SPAdes
        -L  <MIN_LENGTH>      Mimimum contig length.   Default: 300
        -n  <SAMPLE_ID>       Name of the sample.      Default: basename PAIRED_FILE
        "
}

#If less than 2 options are input, show usage and exit script.
if [ $# -le 2 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts s:p:u:l:m:t:e:o:L:n: option ; do
        case "${option}"
        in
                p) PAIRED_FILE=${OPTARG};;
                u) UNPAIRED_FILE=${OPTARG};;
                o) OUTPUT_FILE_NAME=${OPTARG};;

                s) SPADES_TYPE=${OPTARG};;
                l) LOG_FILE=${OPTARG};;
                m) MEMORY=${OPTARG};;
                t) THREADS=${OPTARG};;
                e) TEMP_DIR=${OPTARG};;
                L) MIN_LENGTH=${OPTARG};;
                n) SAMPLE_ID=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Setting defaults
#------------------------------------------------------------------------------#
SPADES_TYPE=${SPADES_TYPE:-rna}
MEMORY=${MEMORY:-10}
THREADS=${THREADS:-1}
TEMP_DIR=${TEMP_DIR:-./SPAdes}
MEGAHIT_ALLOWED=${MEGAHIT_ALLOWED:-no}
MIN_LENGTH=${MIN_LENGTH:-300}
SAMPLE_ID=${SAMPLE_ID:-$(basename $PAIRED_FILE)}

#------------------------------------------------------------------------------#
#Defining functions
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

assembly_with_SPAdes() {
  #Make sure argument values are valid
  if [[ $SPADES_TYPE != "meta" && $SPADES_TYPE != "rna" && $SPADES_TYPE != "regular" ]] ; then
    write_log \
    $SAMPLE_ID \
    "SPADES_TYPE is set to $SPADES_TYPE, but must be 'meta', 'rna', or 'regular'." \
    1 \
    $LOG_FILE
    exit 1
  fi

  #Make sure all directories that are needed are made
  mkdir -p $TEMP_DIR
  mkdir -p $(dirname $OUTPUT_FILE_NAME)

  #Updating log
  write_log \
  $SAMPLE_ID \
  "Starting assembly_with_SPAdes script." \
  0

  #check what combination of files exists, to define which combination
  if [[ -s $PAIRED_FILE ]] && [[ -s $UNPAIRED_FILE ]] ; then
    local PAIRING=BOTH
    write_log $SAMPLE_ID "Found both a paired and unpaied file." 0
  elif [[ -s $PAIRED_FILE ]] && [[ ! -s $UNPAIRED_FILE ]] ; then
    local PAIRING=PAIRED
    write_log $SAMPLE_ID "Found Only a paired file." 0
  elif [[ ! -s $PAIRED_FILE ]] && [[ -s $UNPAIRED_FILE ]] ; then
    local PAIRING=UNPAIRED write_log $SAMPLE_ID "Found only an unpaied file." 0
  elif [[ ! -s $PAIRED_FILE ]] && [[ ! -s $UNPAIRED_FILE ]] ; then
    write_log $SAMPLE_ID "Cannot find either a paired or unpaired file!" 1 $LOG_FILE
    exit 1
  fi

  #defining the SPAdes type
  if [[ $SPADES_TYPE == "meta" ]] ; then
    local SPADES=metaspades.py
  elif [[ $SPADES_TYPE == "rna" ]] ; then
    local SPADES=rnaspades.py
  elif [[ $SPADES_TYPE == "regular" ]] ; then
    local SPADES=spades.py
  fi

  write_log $SAMPLE_ID "SPAdes type is set as $SPADES_TYPE, so running $SPADES." 0

  #Running SPAdes
  if [ $PAIRING == "BOTH" ] ; then
      $SPADES \
      -o $TEMP_DIR/spades \
      --12 $PAIRED_FILE \
      -s $UNPAIRED_FILE \
      -t $THREADS \
      -m $MEMORY \
      --tmp-dir $TEMP_DIR
  elif [ $PAIRING == "PAIRED" ] ; then
      $SPADES \
      -o $TEMP_DIR/spades \
      --12 $PAIRED_FILE \
      -t $THREADS \
      -m $MEMORY \
      --tmp-dir $TEMP_DIR
    elif [ $PAIRING == "UNPAIRED" ] ; then
      $SPADES \
      -o $TEMP_DIR/spades \
      -s $UNPAIRED_FILE \
      -t $THREADS \
      -m $MEMORY \
      --tmp-dir $TEMP_DIR
  fi

  #Determining what output file is expected, and copying it over. Name varies if it's rnaspades.
  if [[ $SPADES_TYPE == "meta" ]] || [[ $SPADES_TYPE == "regular" ]] ; then
    cp $TEMP_DIR/spades/scaffolds.fasta ${OUTPUT_FILE_NAME}.tmp
  elif [[ $SPADES_TYPE == "rna" ]] ; then
    cp $TEMP_DIR/spades/transcripts.fasta ${OUTPUT_FILE_NAME}.tmp
  fi

  # Make sure an outfile was produced by SPAdes. If there is not, this may be
  # a true failure of SPAdes or SPAdes messed up.
  if [[ ! -f ${OUTPUT_FILE_NAME}.tmp ]] ; then
    write_log $SAMPLE_ID "Cannot find the SPAdes output file ${OUTPUT_FILE_NAME}.tmp!\
     Spades may have failed, or this might be a real inability to make contigs." 1 $LOG_FILE
  fi

  # Filter contigs by min-length
  seqtk seq -L $MIN_LENGTH ${OUTPUT_FILE_NAME}.tmp > ${OUTPUT_FILE_NAME} \
    && rm ${OUTPUT_FILE_NAME}.tmp

  # See if the final file has contents. If not, just report it to the log
  if [[ ! -s ${OUTPUT_FILE_NAME} ]] ; then
    write_log $SAMPLE_ID "There are no output contigs following length filtering." 1 $LOG_FILE
  fi

  write_log "Finished assembly_with_SPAdes script for $PAIRED_FILE.:" 0
}


#------------------------------------------------------------------------------#
#Calling functions
#------------------------------------------------------------------------------#
assembly_with_SPAdes
