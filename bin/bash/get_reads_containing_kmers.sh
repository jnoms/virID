#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2


module load gcc

#------------------------------------------------------------------------------#
#Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        This purpose of this file is to find reads containing one or more kmers
        present in the input kmers file, and to isolate those reads to the
        output file.

        Usage:
        $0
        -r <READS_FILE>   Path to the fastq containing all reads.
        -k <KMER_FILE>    Path to file containing kmers, one per one, to search
                          for in the READS_FILE.
        -o <OUTPUT_FILE>  Path to the output fastq containing only reads that
                          have one or more kmers.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts r:k:o: option ; do
        case "${option}"
        in
                r) READS_FILE=${OPTARG};;
                k) KMER_FILE=${OPTARG};;
                o) OUTPUT_FILE=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#

# Search for the kmers in the reads file, and get the headers of positive reads.
# EXPLAINATION - first grep the reads_file for any non-regex (-F) pattern in the
# kmer_file (-f) and report back the line before it (the header). Get only the
# headers by grepping for '@'. Because there are also line numbers, cut those
# out. Use set to remove the "@" from the start of each line.
LC_ALL=C grep -F -b1 -f $KMER_FILE $READS_FILE | grep -F "@" |\
 cut -f1 -d"-" --complement | sed "s/^.//g" > ${OUTPUT_FILE}_headers.txt

 # Pull out identified reads
 seqtk subseq $READS_FILE ${OUTPUT_FILE}_headers.txt > ${OUTPUT_FILE}
 rm ${OUTPUT_FILE}_headers.txt
