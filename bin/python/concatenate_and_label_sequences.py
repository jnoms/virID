#!/usr/bin/env python3

from glob import glob
import sys
import os
import pathlib
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def usage():
    print("""The purpose of this function is to take in a series of either
             fasta's or fastq's, label the reads based on the file names, and
             concatenate them to an output file.

             Usage:
             concatenate_and_label_sequences.py <"input_glob"> <output_path>

             Where,
                 <input_glob> is a glob detailing input files. Inputs must ALL
                 be fasta or ALL be fastq, and must NOT have "__" in the
                 filename. Make sure you're wrapping the glob in quotes.
    """)

    sys.exit()

def get_out_handle(outfile):

    print("Generating output handle")

    #Make directory if necessary
    out_dir = os.path.dirname(outfile)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    #If the file already exists, overwrite it
    if os.path.exists(outfile):
        os.remove(outfile)
    out_handle = open(outfile, "a")

    return out_handle

def validate_infiles(infiles):
    """
    infiles - a list of file pahts

    This function checks that they all end in either fastq or fasta, and
    that each filename does not contain "__", which is what will separate
    file names from sequence headers.
    """

    print("Validating infiles.")

    file_type = ""

    for infile in infiles:

        #Get the file type of the first file
        if file_type == "":
            if infile.endswith("fasta"):
                file_type = "fasta"
                print("Detected fasta files.")

            elif infile.endswith("fastq"):
                file_type = "fastq"
                print("Detected fastq files.")

            else:
                raise ValueError("Infiles must end with .fasta or .fastq")

        #Check that all files are the same type
        if not infile.endswith(file_type):
            raise ValueError("Not all infiles are of the same file type.")

        #Check that the file names don't have a "__" in them
        if "__" in infile:
            raise ValueError("The infile name " + infile +
            " contains a __, which is not allowed.")

def concatenate_fasta(infile, out_handle):
    """
    This function labels each read and add's it to the outfile.

    infile - a simple string with infile path, not a file object
    out_handle - an output file object that has already been "openeded"
    """

    file_name = os.path.basename(infile)[:-len(".fasta")]

    with open(infile) as in_handle:
        for line in in_handle:

            #Add file name to each header
            if line.startswith(">"):
                line = line.replace(">", ">" + file_name + "__")

            #Write the line to the outfile
            out_handle.write(line)

def concatenate_fastq(infile, out_handle):
    """
    This function labels each read and add's it to the outfile.

    infile - a simple string with infile path, not a file object
    out_handle - an output file object that has already been "openeded"
    """

    file_name = os.path.basename(infile)[:-len(".fasta")]

    with open(infile) as in_handle:
        for title, seq, qual in FastqGeneralIterator(in_handle):
            title = file_name + "__" + title
            out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

def main():
    #Check that all args are there
    if len(sys.argv) != 3:
        usage()

    #Prep output file
    out_handle = get_out_handle(sys.argv[2])

    #Prep infiles
    infiles = glob(sys.argv[1])
    validate_infiles(infiles)

    print("Concatenating files.")
    for infile in infiles:
        if infile.endswith("fasta"):
            concatenate_fasta(infile, out_handle)
        elif infile.endswith("fastq"):
            concatenate_fastq(infile, out_handle)

    out_handle.close()

    print("Finished.")

if __name__ == '__main__':
    main()
