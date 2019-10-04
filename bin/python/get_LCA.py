#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 12:29:53 2019

@author: jnom
"""

import sys
minimum_python_version = (3, 6, 7)
if  sys.version_info < minimum_python_version:
    print("""
    This script uses a dictionary I want to be ordered, so it is recommendend
    you use a python version of at least 3.6.7. You can proceed with a lower
    version, but the taxonomy columns may not be in the correct order in the
    output. To proceed anyway, remove this clause from the script.
    """)
    sys.exit()

import inspect
import argparse
import time
import pathlib
import os
from anytree import Node, Walker
import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()
import numpy as np
from collections import Counter

#------------------------------------------------------------------------------#
# Functions
#------------------------------------------------------------------------------#
 # ete3 functions
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

def get_query_ID_counts(infile_contents):
    """ This returns a counter object with the counts for each query_ID."""
    query_ID_list = [line.split('\t')[0] for line in infile_contents]
    counts = Counter(query_ID_list)
    return counts

def read_blacklist_file(blacklist_file_path):
    """
    Input is path to blacklist file. Blacklist file can be either
    structured as <taxonID>\t<any other stuff, i.e. lineage...>
    or just as <taxonID>.

    This function returns a set of taxonIDs that are in the blacklist
    """
    blacklist = set()
    with open(blacklist_file_path, 'r') as infile:
        for line in infile:
            if len(line.split('\t')) >1:
                line = line.split('\t')
                taxonID = line[0]
                blacklist.add(int(taxonID))
            else:
                taxonID = line.rstrip('\n')
                blacklist.add(int(taxonID))

    return blacklist

def taxonID_passes_filtration(line, taxonomy_column_index, current_contents, blacklist):
    """
    Filters based on taxonID. Returns False if the line fails.
    """
    current_taxonIDs = set([item[taxonomy_column_index] for item in current_contents])
    taxonID = line[taxonomy_column_index]

    # taxonID is "N/A", which happens from blast sometimes
    if taxonID == "N/A":
        return False

    # taxonID is blank
    elif taxonID == '':
        return False
        
    # taxonID in blacklist
    elif int(taxonID) in blacklist:
        return False

    # taxonID can't be found in the taxonomy
    elif get_name(taxonID) == "UNKNOWN":
        return False

    # taxonID is already in the current_contents
    elif taxonID in current_taxonIDs:
        return False

    # Otherwise, it passsed
    else:
        return True

def filter_current_contents_by_score(current_contents, score_index, within_percent_of_top_score):

    # Find scores
    scores = [float(item[score_index]) for item in current_contents]

    # Get the threshold, then a boolean for which scores pass the threshold
    max_score = max(scores)
    threshold = (1.0 - (within_percent_of_top_score)/100)*max_score
    score_boolean = np.array([score >= threshold for score in scores])

    # Convert current_contents to a np array so I can do boolean subsetting.
    # Return list.
    current_contents_array = np.array(current_contents)
    return current_contents_array[score_boolean].tolist()

def construct_taxon_tree(taxon_list):
    """
    The purpose of this function is to take an input list of numeric taxonIDs and convert them into a
    dictionary containing anytree nodes with their lineages.

    Input:
    # - taxon_list: A list of numeric taxonIDs

    Output:
    # taxon_tree: A dictionary of node objects, where taxonID:taxon_node_object
    """

    ncbi = NCBITaxa()

    def add_to_taxon_tree(taxon_tree, lineage):
        """
        The purpose of this function is to add a lineage to a taxon tree.

        Input:
        # - taxon_tree: a dictionary containing taxonIDs as keys and taxon anytree node objects as values
        # - lineage: a list of numeric taxonIDs in order of lineage

        Output:
        # - taxon_tree: taxon_tree dictionarywith added lineage if the taxons weren't present in the
            taxon_tree already
        """
        #Initialize tree if needed
        if taxon_tree == dict():
            taxon_tree[1] = Node(1, parent=None)

        for i, taxon in enumerate(lineage):
            #Check if the taxon is already in the tree
            if taxon in taxon_tree:
                continue

            #if taxon is the first in the lineage, parent is 1 (root)
            if i == 0:
                taxon_tree[taxon] = Node(taxon, parent=taxon_tree[1])
                continue

            #otherwise, parent is the previous taxon
            previous_taxon = lineage[i-1]
            taxon_tree[taxon] = Node(taxon, parent=taxon_tree[previous_taxon])

        return taxon_tree

    #State progress
    function_name = inspect.stack()[0][3]
    print(function_name + ': Constructing taxon_tree.')

    taxon_tree = dict()
    for taxon in taxon_list:
        lineage = get_lineage(taxon)

        # Occasionally the terminal taxonID in the lineage isn't the query taxonID (often due to
        # an old reference db). In this case, just reset the terminal taxonID as the query.
        if lineage[-1] != taxon:
            print(function_name + ": The terminal taxonID in the lineage of taxonID " + str(taxon) +
                  " is " + str(lineage[-1]) + ". Overwriting to " + str(taxon))
            lineage[-1] = taxon

        taxon_tree = add_to_taxon_tree(taxon_tree, lineage)

    print(function_name + ': Finished.')
    return taxon_tree

def find_LCA_two_objects(node_object_1, node_object_2):
        w = Walker()
        LCA = w.walk(node_object_1, node_object_2)[1]
        return LCA

def find_LCA_iteratively(taxonIDs, taxon_tree):
    # Initialize LCA
    LCA_lineage = taxon_tree[taxonIDs[0]]

    # find LCA sequentially
    for taxonID in taxonIDs:
        taxon = taxon_tree[taxonID]
        LCA_lineage = find_LCA_two_objects(LCA_lineage, taxon)

    # Format LCA lineage
    LCA_lineage = str(LCA_lineage).rstrip("')'").strip("Node(')")[1:].split('/')

    # Get LCA_taxonID
    LCA_taxonID = int(LCA_lineage[-1])

    return LCA_taxonID


def look_up_LCA_taxonID(taxonIDs):
    """
    This is the parent function. It takes in a list of taxonIDs and returns
    the taxonID of the LCA of those taxonIDs.
    """
    # If taxonIDs is only length of 1, just return the taxonID
    if len(taxonIDs) == 1:
        return taxonIDs[0]

    # Generate a taxon_tree from the input taxonIDs
    taxon_tree = construct_taxon_tree(taxonIDs)

    # Find the taxonID of the LCA
    LCA_taxonID = find_LCA_iteratively(taxonIDs,taxon_tree)

    return LCA_taxonID

def get_cannonical_lineage(taxonID, prefix_dictionary):

    # Get a list of lineage and a list of their ranks
    lineage = get_lineage(taxonID)
    levels = [get_level(taxon) for taxon in lineage]

    cannonical_lineage = []
    # Iterate over each of the levels in prefix dictionary
    for level, prefix in prefix_dictionary.items():


        # If the level isn't here, it is unknown
        if not level in set(levels):
            cannonical_lineage.append(0)
            continue

        # Get the taxon name
        index = levels.index(level)
        taxon = lineage[index]
        name = get_name(taxon)

        # Report it out with the lineage
        cannonical_lineage.append(prefix + name)

    return cannonical_lineage

def aggregate_list(input_list):
    """
    Takes in a list of sublists and reformats. For example

    Input: [[name1, taxonID1, lineage1], [name2, taxonID2, lineage2]]
    Output: [[name1, name2], [taxonID1, taxonID2], [lineage1, lineage2]]
    """
    aggregated = []

    # Iterate as many times as the sublists
    for i in range(len(input_list[0])):
        # We don't care about query_ID bc they're all the same
        if i == 0:
            continue
        aggregated.append([sublist[i] for sublist in input_list])

    return aggregated

def write_output(output, output_file_path):

    #Make output directory if necessary
    output_directory = os.path.dirname(output_file_path)
    pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)

    #Write to output file
    with open(output_file_path, "w") as outfile:
        outfile.write(output)

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

def main():

    #--------------------------------------------------------------------------#
    #Take inputs
    #--------------------------------------------------------------------------#
    parser = argparse.ArgumentParser(description="""
            This program converts a BLAST or DIAMOND output file into report of
            the most likely match.

            If NCBITaxa() from ETE3 has not been
            previously run, ETE3 will download the NCBI taxonomy database to
            ~/.ete3toolkit/. If your home directory is located on a drive with
            slow I/O, I recommend making a symlink from that directory to a
            faster one else speed will suffer. For reference, this script should
            run in about 40s on a BLAST output file of 500K lines.

            Also, to update the taxonomy database run
            ncbi.update_taxonomy_database(). See
            http://etetoolkit.org/docs/latest/tutorial/tutorial_ncbitaxonomy.html
            for more information.
            """)

    # Required arguments
    parser.add_argument(
        '-i',
        '--infile',
        type=str,
        required=True,
        help='''
        Path to the BLAST or DIAMOND output file, that will serve as
        this scripts infile. Assumes the first column contains the sequence ID.
        '''
    )
    parser.add_argument(
        '-o',
        '--outfile',
        type=str,
        required=True,
        help='The path to the output file. File format will be .tsv.'
    )

    #Optional arguments
    parser.add_argument(
        '-c',
        '--column_names',
        type=str,
        required=False,
        default="query_ID seq_title seq_ID taxonID evalue bitscore pident length",
        help='''A space-delimited of column names. Remeber to wrap in quotes.
        <default: query_ID seq_title seq_ID taxonID evalue bitscore pident length>
        '''
    )
    parser.add_argument(
        '-t',
        '--taxonomy_column',
        type=str,
        required=False,
        default="taxonID",
        help='''Name of the column with the taxonIDs.
        <default: taxonID>
        '''
    )
    parser.add_argument(
        '-s',
        '--score_column',
        type=str,
        required=False,
        default="bitscore",
        help='''Name of the column with the scores you want to be filter by.
        <default: bitscore>
        '''
    )
    parser.add_argument(
        '-b',
        '--blacklist_file_path',
        type=str,
        required=False,
        default="",
        help="""
        Path to a file containing taxonIDs which should be blacklisted. Can
        either be a single column with taxonIDs or a tab-delimited file with
        taxonIDs first and anything else afterwards. <default: no blacklist>
        """
    )
    parser.add_argument(
        '-w',
        '--within_percentage_of_top_score',
        type=int,
        required=False,
        default = 1,
        help='''
        Within what percentage of the top score is required to keep. A value of
        1 keeps all hits within 1 percent of top bitscore.
        <default: 1>
        '''
    )
    parser.add_argument(
        '-l',
        '--log_file',
        type=str,
        required=False,
        default = '',
        help='''Path to the log file. If none entered, will print entries to
            stdout. <default: None, prints to screen>''')


    args = parser.parse_args()

    # Define input variables
    infile = args.infile
    outfile = args.outfile
    colnames = args.column_names.split(" ")
    taxonomy_column = args.taxonomy_column
    score_column = args.score_column
    blacklist_file_path = args.blacklist_file_path
    within_percent_of_top_score = args.within_percentage_of_top_score
    log_file = args.log_file

    #--------------------------------------------------------------------------#
    # Constants
    #--------------------------------------------------------------------------#
    prefix_dictionary = {
    'superkingdom': 'sk__',
    'kingdom': 'k__',
    'phylum': 'p__',
    'class': 'c__',
    'order': 'o__',
    'family': 'f__',
    'genus': 'g__',
    'species': 's__',
    'strain': 'st__'
    }
    taxonomy_column_index = colnames.index(taxonomy_column)
    score_index = colnames.index(score_column)

    #--------------------------------------------------------------------------#
    # Main
    #--------------------------------------------------------------------------#

    write_to_log(log_file, "get_LCA.py: Starting.")

    # Read in blacklist if it is specified
    if blacklist_file_path != "":
        blacklist = read_blacklist_file(blacklist_file_path)
    else:
        blacklist = set()

    # Make a string to serve as result. First load up the future column names.
    result = '\t'.join(colnames) + "\tLCA_taxonID\t" + '\t'.join(list(prefix_dictionary.keys())) + '\n'

    # Make holders to hold information until all data for a given query_ID is collected
    current_contents = []

    # Open entire infile to memory
    infile_handle = open(infile)
    infile_contents = infile_handle.read().split('\n')
    infile_handle.close()

    # Keep track of total counts and observed counts
    query_ID_counts = get_query_ID_counts(infile_contents)
    query_ID_observed = dict()

    # Iterate over each line
    for line in infile_contents:

        # Break line into feilds
        line = line.rstrip("\n").split("\t")

        # Get the query_ID
        query_ID = line[0]

        # Add to the count
        query_ID_observed[query_ID] = query_ID_observed.get(query_ID, 0) + 1

        # If line is empty or weird, continue
        if len(line) != len(colnames):
            print("The line " + str(line) + " is not the same length as the colnames, " +
                 " so skipping it.")
            continue

        # Check if taxonID is already in the contents, is blank, or is in the blacklist
        # If it passes, append to current_contents
        if taxonID_passes_filtration(line, taxonomy_column_index, current_contents, blacklist):
            current_contents.append(line)

        # If observed is less than counts, this means there are still more
        # lines with that queryID, so continue.
        if query_ID_observed[query_ID] < query_ID_counts[query_ID]:
            continue

        # If observed is equal to counts, time to do LCA calculation and to export
        elif query_ID_observed[query_ID] == query_ID_counts[query_ID]:

            # If current_contents is empty, that means that none of the lines for that
            # query_ID made it :(. Report that, and continue with loop
            if current_contents == []:
                print("None of the lines from query_ID " + query_ID +
                     " made it through filtration.")
                continue

            # Sanity checks - make sure all the current_contents have the same query_ID
            ID_from_contents = set([item[0] for item in current_contents])
            if len(ID_from_contents) > 1:
                raise ValueError("There is more than one unique query_ID in current contents" +
                                ", which means something is very wrong.")
            if ID_from_contents.pop() != query_ID:
                raise ValueError("The current contents doesn't contain the same query_ID " +
                                "as the current query_ID. Something is wrong.")

            # First, need to filter the current_contents by bitscore
            current_contents = filter_current_contents_by_score(current_contents,
                                                                score_index,
                                                                within_percent_of_top_score)

            # Now, need to get cannonical LCA lineage from a list of taxonIDs
            taxonIDs = [int(item[taxonomy_column_index]) for item in current_contents]
            LCA_taxonID = look_up_LCA_taxonID(taxonIDs)

            # Get the cannonical lineage from the LCA_taxonID
            LCA_cannonical_lineage = get_cannonical_lineage(LCA_taxonID, prefix_dictionary)

            # Aggregate current_contents to group each field in each sublist together
            # Then, make the list of lists a list of comma-delimited strings and
            # replace spaces with underscores
            aggregated = aggregate_list(current_contents)
            aggregated = [",".join(item).replace(" ", "_") for item in aggregated]

            # Add the query_ID, LCA_taxonID and lineage to the aggregated list
            loop_output = [query_ID]
            loop_output.extend(aggregated)
            loop_output.append(LCA_taxonID)
            loop_output.extend(LCA_cannonical_lineage)

            # Make all list items strings, and then condense list into a tab-delimited string
            loop_output = [str(item) for item in loop_output]
            loop_output = "\t".join(loop_output)

            # Make sure the output is the size you'd expect
            if len(loop_output.split('\t')) != len(result.split('\n')[0].split('\t')):
                raise ValueError("The loop output " + str(loop_output) + " is not the same length " +
                                " as the preset colnames in the result variable.")

            # Add this to the result
            result += loop_output + "\n"

            # Wipe current_contents to allow the next iteration of the loop
            current_contents = []

    # Sanity check - make sure the counters are the same
    if query_ID_counts != Counter(query_ID_observed):
        raise ValueError("The query_ID_counts is not the same as the query_ID_observed.")

    write_output(result, outfile)

    write_to_log(log_file, "get_LCA.py: Finished.")


if __name__ == '__main__':
    main()
