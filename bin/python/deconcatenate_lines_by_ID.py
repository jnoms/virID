#!/usr/bin/env python3

import sys
import pathlib

def usage():
    print("""
        The purpose of this function is to read in a tab-delimited stdin,
        isolate the first column, and isolate the sample ID, where
        the first column is <sampleID>__<...some other string...>. It then
        writes each line to a file determined by the sampleID.

        If none of the samples have a sampleID, all sequences are just printed
        to the specified outfile. Please note that if the outfiles exist already
        they will be appended to.

        Usage:
        cat input.txt | deconcatenate_lines_by_ID.py <directory> <out_name>

        Where,
            <directory> Is path to output directory, with or without suffix '/'
            <out_name> Is the default output file name. If there are no
                       sampleID's, output will be <directory>/<out_name>. If
                       there are sampleID's, output will be
                      <directory>/<sampleID>_<out_name>
    """)

    sys.exit()

def remove_sampleID(line):

    #get the first column, and remove the ID which is separated by "__"
    first_col = line.split('\t')[0].split("__")[1]

    #get the remaining columns
    remainder = line.split('\t')[1:]

    #regenerate tab-delimited string
    line = "\t".join([first_col] + remainder)

    return line

def main():

    print('Starting.')

    #Check that all args are there
    if len(sys.argv) != 3:
        usage()

    #Label the args
    directory = sys.argv[1]
    if not directory.endswith("/"):
        directory = directory + "/"
    out_name = sys.argv[2]

    #Make directory if necessary
    pathlib.Path(directory).mkdir(parents=True, exist_ok=True)

    #Determine if there are sampleID's.
    ID_status = ''

    #Making dictionary to store the lines in-memory
    file_contents = dict()

    print('Reading in stdin, sorting into files...')

    #Start reading from stdin
    for line in sys.stdin:

        #Split line into columns
        columns = line.split('\t')

        #Isolate the first column
        first_column = columns[0]

        #Split first column by "__", which should separate sampleID
        header = first_column.split("__")
        if len(header) >2:
            raise ValueError("""There seems to be more than one '__' in the
            fist column.""")

        #Get sampleID
        sampleID = str(header[0])

        #If there is no sampleID...
        if len(header) == 1:

            #Determine if any previous lines had a sampleID
            if ID_status == True:
                raise ValueError("""A previous line had a sampleID but this line
                                    does not. Something is wrong.""")
            ID_status = False

            #Gen output path
            outfile = directory + out_name

            #write to output dict
            file_contents[outfile] = file_contents.get(outfile, '') + line

        #If there is a sampleID...
        if len(header) > 1:

            #Determine if any previous lines did not have a sampleID
            if ID_status == False:
                raise ValueError("""A previous line did not have a sampleID but
                                    this line does. Something is wrong.""")
            ID_status = True

            #Generate output file path
            outfile = directory + sampleID + "_" + out_name

            #Remove the sampleID from the first column if it is there
            line = remove_sampleID(line)

            #write to output dict
            file_contents[outfile] = file_contents.get(outfile, '') + line

    #Write to output files
    print('Writing output files.')
    for sample in file_contents:
        with open(sample, 'a') as outfile:
            outfile.write(file_contents[sample])

    print('Finished.')

if __name__ == '__main__':
    main()
