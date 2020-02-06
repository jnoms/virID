import argparse
import pathlib
import os
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def parse_exclusion_file(exclusion_file, exclusion_column):
    """
    Reads in the specified column of the specified file into a set.
    """

    exclusion_list = set()

    with open(exclusion_file) as infile:
        for line in infile:
            to_exclude = line.split('\t')[exclusion_column]
            exclusion_list.add(to_exclude)

    return exclusion_list

def parse_and_exclude(infile, exclusion_list, fastx_type):
    """
    Takes in path to a fasta or fastq. Keeps sequences whose IDs are not in the
    exclusion_list.
    """

    result = ""

    with open(infile) as infile_handle:

        # Parse based on type
        if fastx_type == "fasta":
            for title, seq in SimpleFastaParser(infile_handle):
                if not title in exclusion_list:
                    result += ">{}\n{}\n".format(title, seq)

        elif fastx_type == "fastq":
            for title, seq, qual in FastqGeneralIterator(infile_handle):
                if not title in exclusion_list:
                    result += "@{}\n{}\n+\n{}\n".format(title, seq, qual)

    return result

def write_output(msg, output_path, append=False):
    """
    Writes the msg to the output path. Makes sure the output_path
    directory exists and all that good stuff.
    """

    # Generate the output directory if necessary
    out_dir = os.path.dirname(output_path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Write output
    if append==True:
        with open(output_path, "a") as outfile_handle:
            outfile_handle.write(msg)
    elif append==False:
        with open(output_path, "w") as outfile_handle:
            outfile_handle.write(msg)


def main():
    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to extract sequences from a fasta or
            fastq that have IDs that are NOT in the specified column of the
            exclusion file. This saves the resultant sequences to the
            designated file. Essentially the opposite of seqtk subseq.
            """)

    # Required arguments
    parser.add_argument(
        '-i',
        '--infile',
        type=str,
        required=True,
        help='''
        Fasta/fastq containing the input sequences.
        '''
    )
    parser.add_argument(
        '-e',
        '--exclusion_file',
        type=str,
        required=True,
        help='''
        Path to exclusion file. Should be tab-delimited. Sequence IDs to NOT
        take from the input_fasta should be in the specified column. Sequnece
        ID's should not have the ">".
        '''
    )
    parser.add_argument(
        '-o',
        '--outfile',
        type=str,
        required=True,
        help='''
        Desired path to the output fasta that lacks the specified reads.
        '''
    )

    # Optional
    parser.add_argument(
        '-f',
        '--fastx_type',
        type=str,
        required=False,
        default="fasta",
        help='''
        Options are "fasta" or "fastq" (default fasta).
        '''
    )
    parser.add_argument(
        '-c',
        '--exclusion_column',
        type=int,
        required=False,
        default=0,
        help='''
        Index of the column that contains the sequence IDs to ignore. This is
        0-indexed, so input accordingly... (default 1)
        '''
    )


    args = parser.parse_args()

    # Define input variables
    infile = args.infile
    exclusion_file = args.exclusion_file
    outfile = args.outfile
    fastx_type = args.fastx_type
    exclusion_column = args.exclusion_column

    # Validate input
    if not fastx_type in {"fastq", "fasta"}:
        msg = "fastx_type must be 'fastq' or 'fasta'. You entered {}.".format(fastx_type)
        raise ValueError(msg)

    #--------------------------------------------------------------------------#
    # Main
    #--------------------------------------------------------------------------#

    # First, read in exclusion file and make exclusion list
    print("Starting. Reading in {}".format(exclusion_file))
    exclusion_list = parse_exclusion_file(exclusion_file, exclusion_column)

    # Next, read in fasta entries that are in the exclusion_list
    passing_fastx = parse_and_exclude(infile, exclusion_list, fastx_type)

    # Write outputs
    print("Writing out to {}".format(outfile))
    write_output(passing_fastx, outfile, append=False)
    print("Finished.")

if __name__ == '__main__':
    main()
