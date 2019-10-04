#!/usr/bin/env python3

import os
import argparse
import pathlib
import pandas as pd

def main():
    #--------------------------------------------------------------------------#
    # Parse Arguments
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
    The purpose of this script is to divide the 'count' column by a specified
    divisor, and multiply by a specified multiplier. Thus,

    New_col_value = (old_col_value * multiplier)/divisor
    """)

    parser.add_argument(
        '-i',
        '--infile',
        type=str,
        required=True,
        help="""
        The input table. This script assumes it is tab-delimieted.
        """
    )
    parser.add_argument(
        '-o',
        '--outfile',
        type=str,
        required=True,
        help="""
        Path to desired output file. Any required dirs will be made.
        """
    )

    parser.add_argument(
        '-ho',
        '--host_reads',
        type=int,
        required=True,
        help="""
        Interger value you want to divide the specified column by.
        """
    )

    parser.add_argument(
        '-m',
        '--multiplier',
        type=int,
        required=False,
        default=1000000,
        help="""
        Value you want to multiply by. Default is 1M. I.e. if you want
        reads per 1M host reads, and the divisor is number of host reads, this
        param should be 1M.
        """
    )

    args = parser.parse_args()
    infile_path = args.infile
    outfile_path = args.outfile
    host_reads = args.host_reads
    multiplier = args.multiplier

    #--------------------------------------------------------------------------#
    # Main
    #--------------------------------------------------------------------------#

    # Open up infile
    infile = pd.read_csv(infile_path, sep='\t', index_col='taxonID')

    # Modify count column
    infile['count'] = infile['count']*multiplier/host_reads

    # Make the output directory if necessary
    output_directory = os.path.dirname(outfile_path)
    pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)

    # Write to outfile
    infile.to_csv(outfile_path, sep='\t')

if __name__ == '__main__':
    main()
