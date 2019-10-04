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
                               something like SAMPLE'_congtigs.fasta.'

        Optional:
        -s  <SPADES_TYPE>     meta, rna, or regular.   Default: rna
        -l  <LOG_FILE>        Log file.                Default: no log file
        -m  <MEMORY>          Number of GB of memory.  Default: 10
        -t  <THREADS>         Number of threads        Default: 1
        -e  <TEMP_DIR>        TEMP directory.          Default: ./SPAdes
        -g  <MEGAHIT_ALLOWED> Allow MEGAHIT to run for samples that fail to
                              assemble through SPAdes? Options are 'yes' or 'no'
                                                       Default: no
        "
}

#If less than 2 options are input, show usage and exit script.
if [ $# -le 2 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts s:p:u:l:m:t:e:o:g: option ; do
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
                g) MEGAHIT_ALLOWED=${OPTARG};;
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

assembly_with_SPAdes() {
  #Make sure argument values are valid
  if [[ $SPADES_TYPE != "meta" && $SPADES_TYPE != "rna" && $SPADES_TYPE != "regular" ]] ; then
    write_log "assembly_with_SPAdes.bash, $SAMPLE: SPADES_TYPE is set to $SPADES_TYPE, but must be 'meta', 'rna', or 'regular'." $LOG_FILE "error"
    exit 1
  fi

  #Make sure all directories that are needed are made
  mkdir -p $(dirname $LOG_FILE )
  mkdir -p $TEMP_DIR
  mkdir -p $(dirname $OUTPUT_FILE_NAME)

  #Get sample name
  local SAMPLE=$( basename $PAIRED_FILE )

  #Updating log
  write_log "assembly_with_SPAdes.bash, $SAMPLE: Starting assembly_with_SPAdes script." $LOG_FILE

  #check what combination of files exists, to define which combination
  if [[ -s $PAIRED_FILE ]] && [[ -s $UNPAIRED_FILE ]] ; then
    local PAIRING=BOTH
    write_log "assembly_with_SPAdes.bash, $SAMPLE: Found both a paired and unpaired file." $LOG_FILE
  elif [[ -s $PAIRED_FILE ]] && [[ ! -s $UNPAIRED_FILE ]] ; then
    local PAIRING=PAIRED
    write_log "assembly_with_SPAdes.bash, $SAMPLE: Found only a PAIRED file." $LOG_FILE
  elif [[ ! -s $PAIRED_FILE ]] && [[ -s $UNPAIRED_FILE ]] ; then
    local PAIRING=UNPAIRED
    write_log "assembly_with_SPAdes.bash, $SAMPLE: Found only an UPAIRED file." $LOG_FILE
  elif [[ ! -s $PAIRED_FILE ]] && [[ ! -s $UNPAIRED_FILE ]] ; then
    write_log "assembly_with_SPAdes.bash, $SAMPLE, ERROR: Cannot find PAIRED or UNPAIRED file. This should be impossible." $LOG_FILE "error"
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
  write_log "assembly_with_SPAdes.bash, $SAMPLE: SPAdes type is set as $SPADES_TYPE, so running $SPADES" $LOG_FILE

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
    cp $TEMP_DIR/spades/scaffolds.fasta $OUTPUT_FILE_NAME
  elif [[ $SPADES_TYPE == "rna" ]] ; then
    cp $TEMP_DIR/spades/transcripts.fasta $OUTPUT_FILE_NAME
  fi

  #Making sure there is an output file. If not, run MEGAHIT if allowed. If not allowed, then exit script.
  if [[ ! -s $OUTPUT_FILE_NAME ]] && [[ $MEGAHIT_ALLOWED == "yes" ]] ; then
      write_log "assembly_with_SPAdes.bash, $SAMPLE, ERROR: Cannot find a spades output file for $PAIRED_FILE. Make sure to check on this file. Running MEGAHIT." $LOG_FILE "error"
      RUN_MEGAHIT
  elif [[ ! -s $OUTPUT_FILE_NAME ]] && [[ $MEGAHIT_ALLOWED == "no" ]] ; then
      write_log "assembly_with_SPAdes.bash, $SAMPLE, ERROR: Cannot find a spades output file for $PAIRED_FILE. Make sure to check on this file. Exiting." $LOG_FILE "error"
      exit 1
  fi

  write_log "assembly_with_SPAdes.bash, $SAMPLE: Finished assembly_with_SPAdes script for $PAIRED_FILE. Time is:" $LOG_FILE
}

RUN_MEGAHIT() {
  #########################
  #DESCRIPTION:
  #########################
  #This function serves to run MEGAHIT assembly.
  #Required inputs are PAIRED_FILE, UNPAIRED_FILE, TEMP_DIR, THREADS, LOG_FILE, and OUTPUT_FILE_NAME.
  #MEGAHIT automatically uses 90% of available memory.

  #########################
  #Starting script:
  #########################
  #Make sure the necessary directories have been made
  mkdir -p $TEMP_DIR
  mkdir -p $(dirname $OUTPUT_FILE_NAME)
  mkdir -p $(dirname $LOG_FILE)

  #Defining SAMPLE
  local SAMPLE=$( basename $PAIRED_FILE )

  #Reporting start of the function
  write_log "assembly_with_SPAdes.bash, $SAMPLE, xxxxMEGAHITxxxx: Starting megahit." $LOG_FILE

  #Determine if the input is paired or unpaired
  if [[ -s $PAIRED_FILE ]] && [[ -s $UNPAIRED_FILE ]] ; then
    local PAIRING=BOTH
    write_log "assembly_with_SPAdes.bash, $SAMPLE, xxxxMEGAHITxxxx: Found both a paired and unpaired file." $LOG_FILE
  elif [[ -s $PAIRED_FILE ]] && [[ ! -s $UNPAIRED_FILE ]] ; then
    local PAIRING=PAIRED
    write_log "assembly_with_SPAdes.bash, $SAMPLE, xxxxMEGAHITxxxx: Found only a paired file." $LOG_FILE
  elif [[ ! -s $PAIRED_FILE ]] && [[ -s $UNPAIRED_FILE ]] ; then
    local PAIRING=UNPAIRED
    write_log "assembly_with_SPAdes.bash, $SAMPLE, xxxxMEGAHITxxxx: Found only an upaired file." $LOG_FILE
  elif [[ ! -s $PAIRED_FILE ]] && [[ ! -s $UNPAIRED_FILE ]] ; then
    write_log "assembly_with_SPAdes.bash, $SAMPLE, xxxxMEGAHITxxxx: Cannot find a paired or unpaired file. Exiting." $LOG_FILE "error"
    exit 1
  fi

  #Running MEGAHIT. Will output contigs to temporary directory, and then will copy over to final directory.
  echo "assembly_with_SPAdes.bash, $SAMPLE, xxxxMEGAHITxxxx: Running MEGAHIT." >> $LOG_FILE
  if [ $PAIRING == "BOTH" ] ; then
    megahit \
    --12 $PAIRED_FILE \
    -r $UNPAIRED_FILE \
    -t $THREADS \
    --presets meta-sensitive \
    -o $TEMP_DIR/MEGAHIT \
    --out-prefix $SAMPLE'_MEGAHIT_OUTPUT'
  elif [ $PAIRING == "PAIRED" ] ; then
    megahit \
    --12 $PAIRED_FILE \
    -t $THREADS \
    --presets meta-sensitive \
    -o $TEMP_DIR/MEGAHIT \
    --out-prefix $SAMPLE'_MEGAHIT_OUTPUT'
  elif [ $PAIRING == "UNPAIRED" ] ; then
    megahit \
    -r $UNPAIRED_FILE \
    -t $THREADS \
    --presets meta-sensitive \
    -o $TEMP_DIR/MEGAHIT \
    --out-prefix $SAMPLE'_MEGAHIT_OUTPUT'
  fi

  #Copying output file to final location. Removing the spaces and making them underscores.
  rm $OUTPUT_FILE_NAME
  cat $TEMP_DIR/MEGAHIT/$SAMPLE'_MEGAHIT_OUTPUT.contigs.fa' | while read LINE; do LINE=${LINE// /_} ; echo $LINE >> $OUTPUT_FILE_NAME ; done

  #Make sure the outfile exists and has contents
  if [[ ! -s $OUTPUT_FILE_NAME ]] ; then
    write_log "assembly_with_SPAdes.bash, $SAMPLE, ERROR: xxxxMEGAHITxxxx: Cannot find a spades output file for $PAIRED_FILE, paired and/or unpaired. Make sure to check on this file. Exiting." $LOG_FILE "error"
    exit 1
  fi

  #Reporting end of the function
  write_log "assembly_with_SPAdes.bash, $SAMPLE, xxxxMEGAHITxxxx: Finished megahit. Time is " $LOG_FILE
}

#------------------------------------------------------------------------------#
#Calling functions
#------------------------------------------------------------------------------#
assembly_with_SPAdes
