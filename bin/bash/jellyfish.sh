#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2


#------------------------------------------------------------------------------#
#Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        This purpose of this file is to use jellyfish2 to extract kmers of a
        specified length from an input fastq.

        Usage:
        $0

        Required:
        -i <INPUT_FILE>   Path to input file
        -o <OUTPUT_FILE>  Path to output file (fasta format).
        -k <KMER_SIZE>    Size of kmers to extract. Consider 21.
        -t <THREADS>      Number of threads to use

        Optional:
        -L <LOW_COUNT> Lowest count kmer that is output. Default is 2, which
                       means no kmers are outputed that are lower frequency than
                       2.
        "
}

#If less than 4 options are input, show usage and exit script.
if [ $# -le 4 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts i:o:k:t:L: option ; do
        case "${option}"
        in
                i) INPUT_FILE=${OPTARG};;
                o) OUTPUT_FILE=${OPTARG};;
                k) KMER_SIZE=${OPTARG};;
                t) THREADS=${OPTARG};;
                L) LOW_COUNT=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#

# Load required dependencies
module load gcc jellyfish

# Define defaults
LOW_COUNT=${LOW_COUNT:-2}

# Run jellyfish
jellyfish count -m $KMER_SIZE -s 100M -t $THREADS -C $INPUT_FILE -o ${OUTPUT_FILE}.jf

# Dump kmers to fasta
jellyfish dump -L $LOW_COUNT ${OUTPUT_FILE}.jf > ${OUTPUT_FILE}
