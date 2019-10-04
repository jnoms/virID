#!/usr/bin/env python3

from glob import glob
from Bio.Seq import Seq
import os
import argparse
import time
import pathlib

#------------------------------------------------------------------------------#
# Define functions
#------------------------------------------------------------------------------#
def write_to_log(log_file, message):

    current_date_time = time.strftime("%c")

    #if the log_file wasn't input, just print message to screen
    if log_file == '':
        print(message)
        print(current_date_time)
        return None

    #Check if the file exists. If not, make the directory if necessary
    log_file_directory = os.path.dirname(log_file)
    pathlib.Path(log_file_directory).mkdir(parents=True, exist_ok=True)

    #Open the log_file in append mode and write message to it
    with open(log_file, 'a') as infile:
        infile.write(message + '\n')
        infile.write(current_date_time + '\n')

def extract_kmers(infile_path):
    output_set = set()
    with open(infile_path) as infile:
        for line in infile:
            if line.startswith(">"):
                continue

            kmer = line.rstrip('\n')

            output_set.add(kmer)

            # also add it's rev compliment
            rev_compliment = str(Seq(kmer).reverse_complement())
            output_set.add(rev_compliment)

    return output_set

def get_dict_from_infile_list(infile_list):
    """
    infile_list is a list of infiles from which to extract kmers and add to the
    output dictionary.
    """
    out_dict = dict()
    out_set = set()
    for infile_path in infile_list:
        print('Processing the infile ' + infile_path)
        sample_name = os.path.basename(infile_path).split(".")[0]
        kmer_set = extract_kmers(infile_path)
        out_dict[sample_name] = kmer_set
        out_set.update(kmer_set)
    return out_dict, out_set
#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
def main():
    parser = argparse.ArgumentParser(description="""
    The purpose of this script is to find which kmers are enriched in
    experimental samples compared to control samples.

    This script takes in globs specifying fastas that contain kmers present in
    experimentals and controls. It then identifies kmers that are present in at
    least EXPERIMENTAL_THRESHOLD of the experimental files, and determines which
    of those kmers are not present in any of the control samples.
    """)

    parser.add_argument(
        '-e',
        '--experimental_files_glob',
        type=str,
        required=True,
        help="""
        A glob that specifies fastas containing the kmers from experimental
        files.
        """
    )
    parser.add_argument(
        '-c',
        '--control_files_glob',
        type=str,
        required=True,
        help="""
        A glob that specifies fastas containing the kmers from control
        files.
        """
    )
    parser.add_argument(
        '-o',
        '--outfile_path',
        type=str,
        required=True,
        help="""
        Path to the output file, which is a txt file where each line is one
        enriched kmer.
        """
    )

    parser.add_argument(
        '-t',
        '--threshold',
        type=int,
        required=True,
        help="""
        The number of experimental samples that should have the kmer for it to
        be considered for enrichment.
        """
        )

    parser.add_argument(
        '-l',
        '--log_file',
        type=str,
        required=False,
        default="",
        help="""
        Path to the log file. If blank, log will be written to screen. Log is in
        tsv format.
        """
        )

    args = parser.parse_args()
    experimental_files_glob = args.experimental_files_glob
    control_files_glob = args.control_files_glob
    outfile_path = args.outfile_path
    threshold = args.threshold
    log_file = args.log_file

    # Gen dictionaries, where sampleID:kmers, as well as all the kmers in a set
    experimental_dict, experimental_kmers = get_dict_from_infile_list(glob(experimental_files_glob))
    control_dict, control_kmers = get_dict_from_infile_list(glob(control_files_glob))
    write_to_log(log_file, "experimental_group_total_kmers\t" + str(len(experimental_kmers)))
    write_to_log(log_file, "control_group_total_kmers\t" + str(len(control_kmers)))

    # For each of the kmers, check how many samples' the kmer is present in and
    # get a count. If it is above threshold, add to common_experimental_kmers.
    common_experimental_kmers = set()
    for kmer in experimental_kmers:
        count = 0
        for sample, kmers in experimental_dict.items():
            if count >= threshold:
                continue
            if kmer in kmers:
                count += 1
        if count >= threshold:
            common_experimental_kmers.add(kmer)
    write_to_log(log_file, "common_experimental_kmers\t" + str(len(common_experimental_kmers)))

    # Determine which of the common experimental kmers are not present in the control kmers
    enriched_kmers = common_experimental_kmers - control_kmers
    write_to_log(log_file, "enriched_kmers\t" + str(len(enriched_kmers)))

    # write enriched kmers to outfile
    with open(outfile_path, "w") as outfile:
        for kmer in enriched_kmers:
            outfile.write(kmer + '\n')

if __name__ == '__main__':
    main()
