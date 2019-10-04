#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 18:40:13 2019

@author: jnom
"""

import pandas as pd
import inspect
import argparse
import time
import os
import pathlib
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import numpy as np


"""
This script has functions to do the following

    1) Add unmarked contig-names to the BLASTN/DIAMOND
       LCA outputs:
           get_headers_from_fasta()
           mark_unassigned_sequences()

    2) Add counts information to the BLASTN/DIAMOND files:
           make_counts_dictinoary_from_counts_file()
           add_read_counts_to_dataframe

    ** The functions above can all be called with
       add_unassigned_and_counts_to_dataframe()


    3) Merge DIAMOND/BLASTN LCA dataframes
            merge_dataframes()

"""
#------------------------------------------------------------------------------#
# Defining Functions
#------------------------------------------------------------------------------#
#Ete3 functions
def get_level(taxonID):
    level = list(ncbi.get_rank([taxonID]).values())

    #Unknown taxonID would yield [], which can't be indexed by [0] to get the string
    if level == []:
        level = "UNKNOWN"
    else:
        level = level[0]
    return level

def get_name(taxonID):
    name = list(ncbi.get_taxid_translator([taxonID]).values())

    #Unknown taxonID would yield [], which can't be indexed by [0] to get the string
    if name == []:
        name = "UNKNOWN"
    else:
        name = name[0]

    name = name.replace(" ", "_")
    return name

def get_lineage(taxonID):
    try:
        lineage = ncbi.get_lineage(taxonID)
    except ValueError:
        print("Cannot find taxonID " + str(taxonID))
        lineage = [taxonID]
    if lineage == None:
        lineage = [0]

    return lineage

def get_named_lineage_from_lineage(lineage, prefix_dictionary =
                                {'strain': 'st__',
                                'species': 's__',
                                'genus': 'g__',
                                'family': 'f__',
                                'order': 'o__',
                                'class': 'c__',
                                'phylum': 'p__',
                                'kingdom': 'k__',
                                'superkingdom': 'sk__'}
):

    named_lineage = []
    conversion_dict = ncbi.get_taxid_translator(lineage)

    for taxonID in lineage:
        level = get_level(taxonID)
        prefix = prefix_dictionary.get(level, 'UNKNOWN__')
        name = prefix + conversion_dict.get(taxonID, "UNKNOWN").replace(" ", "_")
        named_lineage.append(name)

    return named_lineage


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

def read_data_file(input_data_file, sep = '\t', column_names = ''):
    """
    This function reads in a tab-delimited file and assigns the column names appropriately.
    If column names are provided they're used, else it uses the first line of the file as the header.
    column_names is a string that is a space-delimited list.
    """

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Reading the input file ' + input_data_file)

    #Main
    if column_names == '':
        DF = pd.read_csv(input_data_file, engine = 'python', sep = sep, header = 0).fillna(0)
    else:
        DF = pd.read_csv(input_data_file, engine = 'python', sep = sep, header = None).fillna(0)
        DF.columns = column_names.split(" ")

    #Make sure the LCA_taxonID column is the correct datatype
    DF['LCA_taxonID'] = DF['LCA_taxonID'].astype(int)

    print(function_name + ': Finished.')

    return DF

def get_headers_from_fasta(input_fasta):
    '''
    The purpose of this function is to read a fasta and extract an ordered list of fasta headers.

    #Input:
    # - <input_fasta> - Path to the fasta file to be read

    #Output:
    # - <headers> - A list of read headers
    '''

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ": Reading fasta headers from " + input_fasta)

    #Main
    headers = []
    with open(input_fasta, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                header = line.rstrip('\n')[1:]
                headers.append(header)

    #Make sure the headers are all unique
    if len(headers) != len(set(headers)):
        raise ValueError(function_name + ": The headers of the fasta are not all unique.")

    print(function_name + ': Finished.')

    return headers

def mark_unassigned_sequences(input_DF, headers):
    '''
    The purpose of this function is to determine which contigs are missing in the input_DF. input_DF must have a column
    called 'query_ID' and a column named 'superkingdom'. Headers is a list of fasta read headers.
    '''

    def column_names_in_DF(input_DF, required_column_names):
        '''
        The purpose of this function is to determine if each column name in a list of required_column_names
        is present in the input_DF.

        #Input:
        # - input_DF - an input DF
        # - required_column_names - a list of column names that are required

        #Output:
        # - True if the required column names are in the input_DF, false if not.
        '''
        output = True
        for column_name in required_column_names:
            if column_name not in input_DF.columns:
                output = False
        return output

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Writing a new dataframe that contains contigs that were not assigned.')

    #Make sure input_DF has a column named 'query_ID' and 'superkingdom'
    if not column_names_in_DF(input_DF, ['query_ID', 'superkingdom']):
        raise ValueError(function_name + ": The column names 'query_ID' and 'superkingdom' are required.")

    #Find a list of query_IDs present in the input_DF, and make sure they are unique.
    input_DF_query_IDs = list(input_DF['query_ID'].unique())
    if input_DF_query_IDs != list(input_DF['query_ID']):
        raise ValueError(function_name + ": The query_ID's present in the input_DF are not unique.")

    #Make an output DF
    output_DF = pd.DataFrame(columns = input_DF.columns)

    #Iterate over the query_IDs in the headers list and, if the query_ID is in the input_DF, take it.
    #If not, make an empty entry.
    for query_ID in headers:

        #Check if the contig is in the input dataframe
        if query_ID in input_DF_query_IDs:
            output_DF = output_DF.append(input_DF[input_DF['query_ID'] == query_ID])

        else:
            output_DF = output_DF.append(pd.Series(), ignore_index = True)
            output_DF.iloc[-1, output_DF.columns.get_loc('query_ID')] = query_ID
            output_DF.iloc[-1, output_DF.columns.get_loc('superkingdom')] = 'UNASSIGNED'

    if len(output_DF) != len(headers):
        raise ValueError(function_name + ": The length of the output dataframe is not equal to the length of the input contig name list. Something is wrong.")

    print(function_name + ': Finished.')
    return output_DF

def make_counts_dictionary_from_counts_file(counts_file, sep = '\t'):
    """
    The purpose of this function is to read in a counts file and make it a dictionary.
    The structure of the counts file must be a two column file with
    'query_ID' and 'counts' as the columns. There should be no colnames in the file.
    """

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Reading in the counts file to make a counts dictionary.')

    counts_dictionary = dict()
    with open(counts_file, 'r') as infile:
        for line in infile:
            line = line.split(sep)

            if len(line) != 2:
                raise ValueError(function_name + ": The counts_file should have two columns. Make sure the delimiter, sep, is correct.")

            query_ID = line[0]
            count = int(line[1])

            if query_ID == '*':
                continue

            counts_dictionary[query_ID] = count

    print(function_name + ': Finished.')
    return counts_dictionary

def add_read_counts_to_dataframe(input_DF, counts_dictionary):

    '''
    The purpose of this function is to add a column, labeled counts, to the input dataframe. This assumes that
    the dictionary keys of the counts_dictionary are equal to values in the query_ID column of the input dataframe.
    '''

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Adding read counts to dataframe.')

    #Make an output DF that is a copy of the input dataframe
    output_DF = input_DF.copy()

    #Get a list of query_IDs in the input_DF and make sure they are unique
    output_DF_query_IDs = list(output_DF['query_ID'].unique())
    if output_DF_query_IDs != list(output_DF['query_ID']):
        raise ValueError(function_name + ": The query_ID's present in the input_DF are not unique.")

    # If counts dictionary is empty, just set read_counts to 1 (assuming these are reads)
    if counts_dictionary == '':
        output_DF['read_count'] = 1
        print(function_name + ' Did not detect a counts dictionary. Assuming this is reads.')
        print(function_name + ': Finished.')
        return output_DF

    #Make sure there is a count for every entry in the input_DF
    if set(output_DF_query_IDs) != set(counts_dictionary.keys()):
        raise ValueError(function_name + "The contig names in the dataframe and those in the counts dictionary are not the same.")

    #Add a read_count column to the input dataframe
    output_DF['read_count'] = 0

    #For each query_ID, add counts to the output dictionary
    for query_ID in output_DF_query_IDs:
        output_DF['read_count'][output_DF['query_ID'] == query_ID] = counts_dictionary[query_ID]

    print(function_name + ': Finished.')
    return output_DF

def add_unassigned_and_counts_to_dataframe(contig_taxonomy_file, contig_fasta, counts_file):
    """
    The point of this function is to stitch together functions to add unassigned contigs, as well as contig counts,
    to the input contig_taxonomy_file.

    # Function dependencies:
    # - read_data_file
    # - get_headers_from_fasta
    # - mark_unassigned_sequences
    # - make_counts_dictionary_from_counts_file
    # - add_read_counts_to_dataframe

    # INPUT:
    # - [contig_taxonomy_file] - Path to the contig_taxonomy_file for a given sample. This file must have a column named
        query_ID and one named superkingdom. It is expected that the column names are in the first row.
    # - [contig_fasta] - Path to the fasta containig all contigs.
    # - [counts_file] - Path to the file containing counts for each contig. The structure of the counts file must be a two column file with
        'query_ID' and 'counts' as the columns. There should be no colnames in the file.

    # OUTPUT:
    # - output_DF - Contains a pandas dataframe of the input contig_taxonomy_file, but marks which contigs were
        unassigned and gives each contig read counts
    """

    #read in contig_taxonomy file
    raw_data = read_data_file(contig_taxonomy_file)

    #get headers from fasta
    headers = get_headers_from_fasta(contig_fasta)

    #fill in the unassigned contigs
    data = mark_unassigned_sequences(raw_data, headers)

    #Make counts_dictionary from counts_file
    if counts_file == '':
        counts_dictionary = ''
        print("Did not detect a counts file! Assuming these are reads.")
    else:
        counts_dictionary = make_counts_dictionary_from_counts_file(counts_file)

    #Add counts to the dataframe
    data_with_counts = add_read_counts_to_dataframe(data, counts_dictionary)

    return data_with_counts


def merge_dataframes(DF1, DF2, origins):
    """
    The purpose of this function is to take two DF's, merge them based on contents of the query_ID column.

    # Input:
    # - DF1, DF2 - Pandas dataframes containing contigs. 'query_ID' column is required. The contig names present
                   in this column must be the same for each of the two dataframes
    # - origins - A list of length two containing origins in the following structure - ['DF1_origin', 'DF2_origin'].
                  An example input is ['BLASTN', 'DIAMOND'] if DF1 is from BLASTN and DF2 is from DIAMOND.
    # - colnames - A list of columns that are desired to be first in the output dataframe. Should be a python list.
                        A default is set for a common use-case.

    # Output:
    # - A dataframe containing the contents of both input dataframes, ordered by query_ID.
    """

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Merging dataframes.')

    #Make copies of input DF's so I'm not mutating the inputs
    DF1 = DF1.copy()
    DF2 = DF2.copy()

    #Make sure the origins list is of length two
    if len(origins) != 2:
        raise ValueError(function_name + ": origins variable must be a list of two objects.")

    #Add origin to each dataframe
    DF1['origin'] = origins[0]
    DF2['origin'] = origins[1]

    #Make sure the input DF's are of the same length
    if len(DF1) != len(DF2):
        raise ValueError(function_name + ": The input dataframes are not of the same length. Are all contigs \
                         present in both?")

    #make sure both input DF's have a query_ID column
    if 'query_ID' not in (list(DF1.columns) and list(DF2.columns)):
        raise ValueError(function_name + ": query_ID must be a column in both input DFs.")

    # Find colnames. This funny business keeps a sensible colname order.
    total_colnames = list(DF1.columns) + list(DF2.columns)
    new_colnames = []
    for colname in total_colnames:
        if colname not in new_colnames:
            new_colnames.append(colname)

    #Make an empty DF with the colnames
    out_DF = pd.DataFrame(columns = new_colnames)

    #Find a list of query_IDs
    DF1_contigs = list(DF1['query_ID'])
    DF2_contigs = list(DF2['query_ID'])
    if set(DF1_contigs) != set(DF2_contigs):
        raise ValueError(function_name + ": DF1 and DF2 must have the same contig names. Is one DF missing contigs?")

    for contig_list in [DF1_contigs, DF2_contigs]:
        if len(contig_list) != len(set(contig_list)):
            raise ValueError(function_name + ": Not all query_IDs are unique in one of the input DFs.")

    #Iterate through contigs, take the row from each DF and add it to the new DF
    for contig in DF1_contigs:
        out_DF = out_DF.append(DF1[DF1['query_ID'] == contig], sort=True)
        out_DF = out_DF.append(DF2[DF2['query_ID'] == contig], sort=True)

    #Reorder the columns
    out_DF = out_DF[new_colnames]

    print(function_name + ': Finished.')
    return out_DF


def write_output(output_DF, output_path, index=False):
    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Writing output to ' + output_path)

    #Make the output directory if necessary
    output_directory = os.path.dirname(output_path)
    pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)

    #Write the output file
    if index == False:
        output_DF.to_csv(output_path, sep='\t', index=False)
    elif index == True:
        output_DF.to_csv(output_path, sep='\t', index=True)

    print('write_output(): Finished.')


#------------------------------------------------------------------------------#
# Define and execute main
#------------------------------------------------------------------------------#
def main():

    #---------------------------------------------------------------------------#
    # Parse inputs
    #---------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
    The purpose of this script is to add more information to blast and diamond
    taxonomy output files - it marks which sequences were not assigned, and
    writes the read counts to each sequence - and to merge the blast and diamond
    results into a joint dataframe.
    """)

    parser.add_argument(
        "-t",
        "--BLASTN_contig_taxonomy_file",
        type=str,
        required=True,
        help='''Path to the BLASTN get_LCA output file. Column names are read from
        the first row of the file. A column named "query_ID" and one named
        "superkingdom" are required.''',
    )
    parser.add_argument(
        "-T",
        "--DIAMOND_contig_taxonomy_file",
        type=str,
        required=True,
        help='''Path to the DIAMOND get_LCA output file. Column names are read
        from the first row of the file.A column named "query_ID" and one
        named "superkingdom" are required.''',
    )
    parser.add_argument(
        "-f",
        "--contig_fasta_file",
        type=str,
        required=True,
        help="""Path to the .fasta containing the contigs. This is used to see
        which sequences weren't assigned.""",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        required=True,
        help="Path to the merged dataframe output file.",
    )

    # Optional parameters
    parser.add_argument(
        "-c",
        "--contig_counts_file",
        type=str,
        required=False,
        default="",
        help="""Path to the .txt (tab-separated) that contains the counts.
        File is of structure <seq_ID><\t><count>. The usecase if when the
        seq_IDs refer to contigs, and you have the counts for each contig.
        <Default: No counts file. Each sequence assigned count 1>
        """,
    )
    parser.add_argument(
        "-l",
        "--log_file",
        type=str,
        required=False,
        default="",
        help="Path to the log file.",
    )


    args = parser.parse_args()

    BLASTN_contig_taxonomy_file = args.BLASTN_contig_taxonomy_file
    DIAMOND_contig_taxonomy_file = args.DIAMOND_contig_taxonomy_file
    contig_fasta_file = args.contig_fasta_file
    output_file = args.output_file
    contig_counts_file = args.contig_counts_file
    log_file = args.log_file

    #---------------------------------------------------------------------------#
    # MAIN
    #---------------------------------------------------------------------------#
    write_to_log(log_file, "Starting pyscript.")

    BLASTN_df = add_unassigned_and_counts_to_dataframe(BLASTN_contig_taxonomy_file, contig_fasta_file, contig_counts_file)
    DIAMOND_df = add_unassigned_and_counts_to_dataframe(DIAMOND_contig_taxonomy_file, contig_fasta_file, contig_counts_file)
    MERGED_df = merge_dataframes(BLASTN_df, DIAMOND_df, ['BLASTN', 'DIAMOND'])

    #write outputs
    write_output(MERGED_df, output_file)

    write_to_log(log_file, "Finished pyscript.")

if __name__ == '__main__':
    main()
