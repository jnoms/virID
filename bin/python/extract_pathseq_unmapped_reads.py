import pysam
import argparse
import pathlib
import os

def get_unmapped_reads(bam_path):

    bam = pysam.AlignmentFile(bam_path, "rb")
    fastq = ""

    # Iterate through all reads
    for read in bam.fetch(until_eof=True):

        # Keep reads that don't have the YP tag
        if not read.has_tag('YP'):

            read_as_fq = "@{}\n{}\n+\n{}\n".format(read.qname, read.seq, read.qual)
            fastq += read_as_fq

    bam.close()

    return fastq

def write_output(msg, output_path):
    """
    Writes the msg to the output path. Makes sure the output_path
    directory exists and all that good stuff.
    """

    # Generate the output directory if necessary
    out_dir = os.path.dirname(output_path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Write output
    with open(output_path, "w") as outfile_handle:
        outfile_handle.write(msg)

def main():
    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to parse a pathseq output bam and
            convert the reads that dont have the YP tag to fastq - i.e. the
            PathSeq unmapped reads.

            Unfortunately, you can't get PathSeq unmapped reads from the
            PathSeq bam by using samtools and the 4 ("unmapped") flag. This is
            because sometimes a read does not have a 4 flag because it was
            aligned, even though it did not align at a match score above the
            bwa score threshold and was thus not counted in the PathSeq score
            module.

            All reads that were counted in the score have the YP tag.
            """)

    # Required arguments
    parser.add_argument(
        '-b',
        '--bam_path',
        type=str,
        required=True,
        help='''
        Path to the GATK-PathSeq output bam.
        '''
    )
    parser.add_argument(
        '-u',
        '--unmapped_reads_path',
        type=str,
        required=True,
        help='''
        Path to unmapped reads in fastq format.
        '''
    )

    args = parser.parse_args()

    # Define input variables
    bam_path = args.bam_path
    unmapped_reads_path = args.unmapped_reads_path

    #--------------------------------------------------------------------------#
    # Main
    #--------------------------------------------------------------------------#
    print("Parsing {}".format(bam_path))
    unmapped_reads = get_unmapped_reads(bam_path)

    n_unmapped = (len(unmapped_reads.split("\n"))-1)/4
    print("There are {} unmapped reads".format(n_unmapped))

    print("Writing output.")
    write_output(unmapped_reads, unmapped_reads_path)

    print("Finished.")

if __name__ == '__main__':
    main()
