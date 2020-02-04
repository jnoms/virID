import pysam
import argparse
import pathlib
import os

def parse_bamfile(bam_path):
    """
    Given the path to a bamfile, parses the reads and generates
    1) a dictionary of mapped reads, of structure
       reference_name:[(name1, qual1, seq1), (name2, qual2, seq2),...]
    2) A list of unmapped reads, which is just a list of tuples like -
      [(name1, qual1, seq1), (name2, qual2, seq2),...]

    """

    bamfile = pysam.AlignmentFile(bam_path, "r")

    mapped_reads_dict = dict()
    unmapped_reads = []
    readcount = 0
    for line in bamfile:

        # Count number of reads in the bam
        readcount += 1

        # parse bamfile
        read_id = line.qname
        qual = line.qual
        seq = line.seq
        reference_name = line.reference_name

        read_info = (read_id, qual, seq)


        # Write to dictionary or list
        if reference_name != None:

            # First, make entry if it doesn't exist
            if reference_name not in mapped_reads_dict:
                mapped_reads_dict[reference_name] = []

            # Write to mapped reads dict
            mapped_reads_dict[reference_name].append(read_info)

        else:
            unmapped_reads.append(read_info)

    return mapped_reads_dict, unmapped_reads, readcount

def parse_contig_assignments(assignment_file_path, contig_name_col=0, unassigned_col=9):
    """
    Finds contig_names that are unassigned twice.

    Inputs:
    - assignment_file_path: Path to tab-delimited file
    - contig_name_col: Index of the col which will have the contig name
    - unassigned_col: Index of the col which will say if it is unassigned, i.e.
        will say 'UNASSIGNED'

    Output: A set of contig names where there are n unassigned lines.

    """
    traking = set()
    unassigned = set()
    with open(assignment_file_path) as infile_handle:

        for line in infile_handle:
            line = line.split("\t")

            contig_name = line[contig_name_col]
            unassigned_status = line[unassigned_col]

            if unassigned_status == "UNASSIGNED":

                if contig_name in traking:
                    unassigned.add(contig_name)

                else:
                    traking.add(contig_name)

    return unassigned

def move_unassigned_reads_to_unmapped(unassigned_contigs, mapped_reads_dict, unmapped_reads):
    """
    Given the unassigned contigs, moves their reads from mapped_reads_dict and
    appends them to unmapped_reads. Then deletes them from mapped_reads_dict.
    """
    # Make copies of the input
    mapped_reads_dict = mapped_reads_dict.copy()
    unmapped_reads = unmapped_reads.copy()

    for contig in unassigned_contigs:

        # Add to unmapped reads
        unmapped_reads.extend(mapped_reads_dict[contig])

        # Remove from mapped reads dictionary
        del mapped_reads_dict[contig]

    return mapped_reads_dict, unmapped_reads

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

def unmapped_reads_to_fastq_format(unmapped_reads):
    fastq = ""
    for read in unmapped_reads:
        header, qual, seq = read
        read = "@{}\n{}\n+\n{}\n".format(header, seq, qual)
        fastq += read

    return fastq

def main():
    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to read in the virID merged assignment
            file and associated virID mapped bam and 1) write out the
            unassembled reads and the reads that are mapped to unassigned
            contigs to a fasta, and 2) (optional) write out a tab-delimted
            metrics file detailing the total number of reads and the number of
            unmapped+unassigned reads.
            """)

    # Required arguments
    parser.add_argument(
        '-b',
        '--bam_path',
        type=str,
        required=True,
        help='''
        Path to the virID mapped bam, which consists of reads mapped back to the
        contigs.
        '''
    )
    parser.add_argument(
        '-a',
        '--assignment_file_path',
        type=str,
        required=True,
        help='''
        Path to the virID MERGED output assignment file.
        '''
    )
    parser.add_argument(
        '-u',
        '--unassigned_reads_path',
        type=str,
        required=True,
        help='''
        Desired path to the output fastq containing unmapped reads.
        '''
    )

    parser.add_argument(
        '-s',
        '--sample_name',
        type=str,
        required=True,
        help='''
        Name of the sample. Used when writing the metrics file.
        '''
        )

    # Optional
    parser.add_argument(
        '-m',
        '--metrics',
        type=str,
        required=False,
        default="",
        help='''
        Path to the optional metrics file.
        '''
        )



    args = parser.parse_args()

    # Define input variables
    bam_path = args.bam_path
    assignment_file_path = args.assignment_file_path
    unassigned_reads_path = args.unassigned_reads_path
    sample_name = args.sample_name
    metrics = args.metrics

    #--------------------------------------------------------------------------#
    # Main
    #--------------------------------------------------------------------------#
    print("Starting extract_unassigned_reads.py for {}".format(sample_name))

    # First, parse the bam to get a dictionary with mapped reads and a
    # list of unmapped reads.
    mapped_reads_dict, unmapped_reads, input_readcount = parse_bamfile(bam_path)

    # Get contigs which were unassigned by both megablast and diamond
    unassigned_contigs = parse_contig_assignments(assignment_file_path)

    # Move the unassigned reads to the unmapped reads and delete from
    # mapped reads dict
    mapped_reads_dict, unmapped_reads = move_unassigned_reads_to_unmapped(unassigned_contigs,
                                                                          mapped_reads_dict,
                                                                          unmapped_reads)

    # Write unmapped reads to fasta
    unmapped_reads_fasta = unmapped_reads_to_fastq_format(unmapped_reads)
    write_output(unmapped_reads_fasta, unassigned_reads_path, append=False)


    # Write total # and # unmapped reads out to metrics file
    if metrics != "":
        metrics_output = "{}\t{}\t{}\n".format(sample_name, input_readcount, len(unmapped_reads))
        write_output(metrics_output, metrics, append=True)

    print("Finished.")

if __name__ == '__main__':
    main()
