#!/usr/bin/env python3

import pandas as pd
from glob import glob
import pathlib
import inspect
import argparse
import os

#------------------------------------------------------------------------------#
# Defining Constants
#------------------------------------------------------------------------------#
consistent_columns = ['taxonID', 'lineage', 'superkingdom', 'taxon', 'level']

#------------------------------------------------------------------------------#
# Define Functions
#------------------------------------------------------------------------------#
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
        DF = pd.read_csv(input_data_file, engine = 'python', sep = sep, header = 0).fillna(0)
    else:
        DF = pd.read_csv(input_data_file, engine = 'python', sep = sep, header = None).fillna(0)
        DF.columns = column_names.split(" ")

    #Make sure the count column is the correct datatype
    DF['count'] = DF['count'].astype(float)

    print(function_name + ': Finished.')

    return DF


def import_target_dir(input_files_glob):
    """
    Reads in a glob, which should specify multiple files, and returns a list
    of dataframes.
    """
    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ':  Reading in dataframes...')

    df_list = []
    for infile in glob(input_files_glob):

        #get sample name
        sample = pathlib.Path(infile).stem

        #read in input file as a pd dataframe
        df = read_data_file(infile)

        #rename the "count" column to have information on the sample name
        df = df.rename(index=str, columns={"count": sample})


        df_list.append(df)

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ':  Finished. Found ' + str(len(df_list)) +
            " dataframes.")

    return df_list

def aggregate_dataframe_list(dataframe_list):

    def validate_aggregated(aggregated, df_list):
        aggregated_lineages = set(aggregated['lineage'])

        df_list_lineages = set()
        for df in df_list:
            lineages = list(df['lineage'])
            for lineage in lineages:
                df_list_lineages.add(lineage)

        if aggregated_lineages != df_list_lineages:
            raise ValueError("""The aggregated dataframe
            does not have exactly the same total lineages
            as the combined input dataframes. Something
            is wrong.""")

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Aggregating dataframes.')

    def sort_columns(dataframe, first_colnames=[]):
        """
        Sorts columns of dataframe. If you specify first_colnames those will
        go first.
        """
        new_colnames = first_colnames.copy()
        temp_sorted = sorted(list(dataframe.columns))
        for colname in temp_sorted:
            if not colname in new_colnames:
                new_colnames.append(colname)
        dataframe = dataframe.reindex(new_colnames, axis=1)
        return dataframe

    aggregated = ''
    for dataframe in dataframe_list:

        #Load the first dataframe to the aggregated variable if it is not a dataframe
        if not isinstance(aggregated, pd.DataFrame):
            aggregated = dataframe
            continue
        aggregated = pd.merge(aggregated, dataframe, how='outer', on=consistent_columns)

    aggregated = sort_columns(aggregated, ['lineage', 'superkingdom', 'taxon', 'level'])

    validate_aggregated(aggregated, dataframe_list)

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Finished.')

    return aggregated

def write_output(output_DF, output_path):
    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Writing output to ' + output_path)

    #Fill NA's to 0
    output_DF.fillna(value=0, inplace=True)

    #Make the output directory if necessary
    output_directory = os.path.dirname(output_path)
    pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)

    #Write the output file
    output_DF.to_csv(output_path, sep='\t', index=False)

    print('write_output(): Finished.')

#------------------------------------------------------------------------------#
# Define and execute main
#------------------------------------------------------------------------------#
def main():
    parser = argparse.ArgumentParser(description="""
    The purpose of this script is to aggregate input .tsv count files. Input is
    a glob detailing the input files, and path to the desired output file.

    Input:
    # - input_files_glob - A glob detailing input files. Must contain the
        columns 'lineage', 'superkingdom', 'taxon', 'level'. Must be tab-
        delimited.
    # - output_path - A path to the output file.

    Output:
    # - aggregated_dataframe - A dataframe with all read counts. Tab-delimited.
    """)

    parser.add_argument(
        '-i',
        '--input_files_glob',
        type=str,
        required=True,
        help="""
        A glob specifying input files. Must contain the columns 'lineage',
        'superkingdom', 'taxon', 'level'. Must be tab-delimited.
        """
    )
    parser.add_argument(
        '-o',
        '--output_path',
        type=str,
        required=True,
        help="""
        Path to desired output file. Necessary directories will be made if not
        currently extant.
        """
    )

    args = parser.parse_args()

    input_files_glob = args.input_files_glob
    output_path = args.output_path

    #Main
    dataframe_list = import_target_dir(input_files_glob)
    if len(dataframe_list) == 0:
        raise ValueError(
        """
        Failed to find input dataframes. Check that the glob is correct and
        is in quotes.
        """
        )

    aggregated_dataframe = aggregate_dataframe_list(dataframe_list)
    write_output(aggregated_dataframe, output_path)




if __name__ == '__main__':
    main()
