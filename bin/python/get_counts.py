#!/usr/bin/env python3

import pandas as pd
import pathlib
import inspect
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import argparse
import time
import os

#------------------------------------------------------------------------------#
# Defining functions
#------------------------------------------------------------------------------#

#Ete3 functions
def get_level(taxonID):
    level = list(ncbi.get_rank([taxonID]).values())

    # Unknown taxonID would yield [], which can't be indexed by [0] to get the
    # string
    if level == []:
        level = "UNKNOWN"
    else:
        level = level[0]
    return level

def get_name(taxonID):
    name = list(ncbi.get_taxid_translator([taxonID]).values())

    # Unknown taxonID would yield [], which can't be indexed by [0] to get the
    # string
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
    This function reads in a tab-delimited file and assigns the column names
    appropriately. If column names are provided they're used, else it uses the
    first line of the file as the header. column_names is a string that is a
    space-delimited list.
    """

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Reading the input file ' + input_data_file)

    #Main
    if column_names == '':
        DF = pd.read_csv(input_data_file, engine = 'python', sep = sep,
        header = 0).fillna(0)
    else:
        DF = pd.read_csv(input_data_file, engine = 'python', sep = sep,
        header = None).fillna(0)
        DF.columns = column_names.split(" ")

    #Make sure the LCA_taxonID column is the correct datatype
    DF['LCA_taxonID'] = DF['LCA_taxonID'].astype(int)

    print(function_name + ': Finished.')

    return DF

def make_counts_dictionary_from_counts_file(counts_file, sep = '\t'):
    """
    The purpose of this function is to read in a counts file and make it a
    dictionary. The structure of the counts file must be a two column file with
    'query_ID' and 'counts' as the columns. There should be no colnames in
    the file.
    """

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ''': Reading in the counts file to make a counts
        dictionary.''')

    counts_dictionary = dict()
    with open(counts_file, 'r') as infile:
        for line in infile:
            line = line.split(sep)

            if len(line) != 2:
                raise ValueError(function_name + """: The counts_file should
                have two columns. Make sure the delimiter, sep, is correct.""")

            query_ID = line[0]
            count = int(line[1])

            if query_ID == '*':
                continue

            counts_dictionary[query_ID] = count

    print(function_name + ': Finished.')
    return counts_dictionary

def add_read_counts_to_dataframe(input_DF, counts_dictionary):

    '''
    The purpose of this function is to add a column, labeled counts, to the
    input dataframe. This assumes that the dictionary keys of the
    counts_dictionary are equal to values in the query_ID column of the input
    dataframe.
    '''

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Adding read counts to dataframe.')

    #Make an output DF that is a copy of the input dataframe
    output_DF = input_DF.copy()

    #Get a list of query_IDs in the input_DF and make sure they are unique
    output_DF_query_IDs = list(output_DF['query_ID'].unique())
    if output_DF_query_IDs != list(output_DF['query_ID']):
        raise ValueError(function_name + """: The query_ID's present in the
        input_DF are not unique.""")

    #Add a read_count column to the input dataframe
    output_DF['read_count'] = 0

    #For each query_ID, add counts to the output dictionary
    for query_ID in output_DF_query_IDs:
        output_DF['read_count'][output_DF['query_ID'] == query_ID] = counts_dictionary[query_ID]

    print(function_name + ': Finished.')
    return output_DF

def get_taxonID_counts(df):
    '''
    Get total counts for each LCA_taxonID. Required colnames are 'LCA_taxonID'
    and 'read_count'. This takes the sum of read_counts for each taxonID and
    returns a dictionary of structure taxonID:counts. Does NOT do any assigning
    to lineages and whatnot.

    Input:
    # - df: A dataframe that contains, at least, the columns "LCA_taxonID" and
    "read_count"

    Output:
    # - taxon_counts_dict: A dictionary of structure taxonID:read_count, where
     read_count is the sum of all read counts of that taxonID
    '''
    def column_names_in_DF(input_DF, required_column_names):
        output = True
        for column_name in required_column_names:
            if column_name not in input_DF.columns:
                output = False
        return output

    if not column_names_in_DF(df, ['LCA_taxonID', 'read_count']):
        raise ValueError("The columns LCA_taxonID and read_count are required.")

    query_ID_to_taxon = dict(zip(df['query_ID'], df['LCA_taxonID']))
    query_ID_to_counts = dict(zip(df['query_ID'], df['read_count']))

    #Get total counts
    taxon_counts_dict = dict()
    for query_ID in query_ID_to_taxon:
        taxon = query_ID_to_taxon[query_ID]
        counts = query_ID_to_counts[query_ID]
        taxon_counts_dict[taxon] = taxon_counts_dict.get(taxon, 0) + counts

    return taxon_counts_dict

def assign_counts(counts_dict, prefix_dictionary = {'strain': 'st__',
                                                 'species': 's__',
                                                 'genus': 'g__',
                                                 'family': 'f__',
                                                 'order': 'o__',
                                                 'class': 'c__',
                                                 'phylum': 'p__',
                                                 'kingdom': 'k__',
                                                 'superkingdom': 'sk__'}
):

    def distribute_counts(counts_dict):
        read_counts = dict()

        #Iterate over each taxonID
        for LCA_taxonID, counts in counts_dict.items():

            #Check if nan - it will not be equal to itself
            if LCA_taxonID != LCA_taxonID:
                continue

            #Check if taxonID is 0 - no need to count
            if LCA_taxonID == 0:
                continue

            #Get lineage
            lineage = get_lineage(LCA_taxonID)

            #Assign counts to each taxonID in lineage
            for taxonID in lineage:
                read_counts[taxonID] = read_counts.get(taxonID, 0) + counts

        return read_counts

    def generate_output_df(read_counts):
        out_df = pd.DataFrame(columns = ['lineage',
                                        'superkingdom',
                                        'taxon',
                                        'level',
                                        'count'])
        out_df['taxon'] = read_counts.keys()
        out_df['taxon'] = out_df['taxon'].astype(str)
        out_df.index = read_counts.keys()
        out_df.index.name = "taxonID"
        return out_df

    def is_cannonical_level(taxonID, cannonical_levels):
        level = get_level(taxonID)
        if level in cannonical_levels:
            return True
        else:
            return False

    def get_cannonical_lineage(taxonID, cannonical_levels):
        lineage = get_lineage(taxonID)
        cannonical_lineage = []
        for taxon in lineage:
            if is_cannonical_level(taxon, cannonical_levels):
                cannonical_lineage.append(taxon)
        return cannonical_lineage

    def get_superkingdom_from_named_lineage(named_lineage):
        superkingdom = "UNKNOWN"
        for taxon in named_lineage:
            if taxon.startswith("sk__"):
                superkingdom = taxon
                return superkingdom
        return superkingdom

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Generating counts dataframe.')

    cannonical_levels = list(prefix_dictionary.keys())
    read_counts = distribute_counts(counts_dict)

    #generate output_df
    out_df = generate_output_df(read_counts)

    #For each taxon, check if it is a cannonical level. If so, write to output df
    for taxonID, count in read_counts.items():

        #Check if the taxonID is cannonical. If not, continue
        if not is_cannonical_level(taxonID, cannonical_levels):
            continue

        level = get_level(taxonID)
        prefix = prefix_dictionary.get(level, "unknown__")
        name = prefix + get_name(taxonID)
        cannonical_lineage = get_cannonical_lineage(taxonID, cannonical_levels)
        named_lineage = get_named_lineage_from_lineage(cannonical_lineage, prefix_dictionary)
        superkingdom = get_superkingdom_from_named_lineage(named_lineage)

        #Write to output dictionary -
        # 'lineage', 'superkingdom', 'taxon', 'level', 'count'
        out_df.at[taxonID, 'lineage'] = named_lineage
        out_df.at[taxonID, 'superkingdom'] = superkingdom
        out_df.at[taxonID, 'taxon'] = name
        out_df.at[taxonID, 'level'] = level
        out_df.at[taxonID, 'count'] = count

    #Remove rows that have missing data
    out_df = out_df.dropna()

    print(function_name + ': Finished.')

    return out_df

def write_output(output_DF, output_path):
    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Writing output to ' + output_path)

    #Make the output directory if necessary
    output_directory = os.path.dirname(output_path)
    pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)

    #Write the output file
    output_DF.to_csv(output_path, sep='\t', index=True)

    print('write_output(): Finished.')

def main(infile, outfile, counts_file, log_file):

    write_to_log(log_file, "Starting.")

    #import data
    data = read_data_file(infile)

    #if there are no counts, assume these are reads
    if counts_file == '':
        print('''Did not detect a counts file. Assuming the input dataset
              is queried reads.''')
        data['read_count'] = 1
    else:
        #Make counts_dictionary from counts_file
        counts_dictionary = make_counts_dictionary_from_counts_file(counts_file)

        #Add counts to the dataframe
        data = add_read_counts_to_dataframe(data, counts_dictionary)

    #Make a dict containing the aggregate counts for each taxonID
    counts_dict = get_taxonID_counts(data)

    #Assign counts to each taxonID's phylogeny
    counts = assign_counts(counts_dict)

    #write output
    write_output(counts, outfile)

    write_to_log(log_file, "Finished.")

if __name__ == '__main__':

    #Argparse input
    parser = argparse.ArgumentParser(description="""
    The purpose of this script is to take a converted DIAMOND or BLAST output
    file (in csv format) and generate a counts file at each taxon level.

    Input:
    - infile: A converted DIAMOND or blast output file. This script assumes
    there ARE colnames. Required column is 'LCA_taxonID'. If there is no
    counts_file, this input is assumed to be from reads, meaning each line has
    a read_count of 1.
    - outfile: Path to the output file, which will be tab-delimited.

    *Optional*
    - counts_file: A tab-delimited file of structure <sequence name> <count>. If
    this is provided, the counts for each entry in the infile will come from this
    file.
    - log_file: An optional file provided for logging purposes (start/end of the script).

    """)
    parser.add_argument('-i', '--infile', type=str, required=True,
        help='''Path to the input data file. See -h for more information.''')
    parser.add_argument('-o', '--outfile', type=str, required=True,
        help='''Path to the output file (tab-delimited).''')
    parser.add_argument('-c', '--counts_file', type=str, required=False, default = '',
        help='''Path to the counts file. See -h for more information.''')
    parser.add_argument('-l', '--log_file', type=str, required=False, default = '',
        help='''Path to the log file.''')

    args = parser.parse_args()

    infile = args.infile
    outfile = args.outfile
    counts_file = args.counts_file
    log_file = args.log_file

    #Run script.
    main(infile, outfile, counts_file, log_file)
