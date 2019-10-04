#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 12:13:08 2019

@author: jnom
"""

import argparse
import os
import time
import pathlib

print('The purpose of this script is to take read files and split them into  1) paired (interleaved) and 2) unpaired fastqs.')
print('The input for this script is the path to a fastq, and the desired output directory. Use -h for more info.')

parser = argparse.ArgumentParser(description='This program inputs files that may have a mix of paired and unpaired reads, and creates a file containing interleaved paired reads and a file containing unpaired reads.')
parser.add_argument('-f', '--fastq', type=str, required=True, help='Name/path to the input fastq.')
parser.add_argument('-o', '--OUT_PREFIX', type=str, required=True, help='Prefix/path to the output files, where output files will be OUT_PREFIX"_paired.fastq" and OUT_PREFIX"_unpaired.fastq".')
parser.add_argument('-l', '--log', type=str, required=False, default='', help='Path to the log file. Not required.')
args = parser.parse_args()

input_fastq = args.fastq
OUT_PREFIX = args.OUT_PREFIX
log_file = args.log


def write_to_log(log_file, message):

    import time
    import os
    import errno

    current_date_time = time.strftime("%c")

    #if the log_file wasn't input, just print message to screen
    if log_file == '':
        print(message)
        print(current_date_time)
        return None


    #Make the log file directory if necessary
    log_file_directory = os.path.dirname(log_file)
    pathlib.Path(log_file_directory).mkdir(parents=True, exist_ok=True)

    #Open the log_file in append mode and write message to it
    with open(log_file, 'a') as infile:
        infile.write(message + '\n')
        infile.write(current_date_time + '\n')

write_to_log(log_file, "process_read_pairs.py: Starting for " + str(input_fastq))

#Make the output directory if it doesn't exit already
OUT_DIRNAME = os.path.dirname(OUT_PREFIX)
pathlib.Path(OUT_DIRNAME).mkdir(parents=True, exist_ok=True)


print('Detected the infile ' + input_fastq)
print('Prefix is set to ' + OUT_PREFIX)

#Preparing output files
paired_output_file_name = OUT_PREFIX + "_paired.fastq"
paired_output_file = open(paired_output_file_name, "w")
paired_output_file.close()
paired_output_file = open(paired_output_file_name, "a")


unpaired_output_file_name= OUT_PREFIX + "_unpaired.fastq"
unpaired_output_file = open(unpaired_output_file_name, "w")
unpaired_output_file.close()
unpaired_output_file = open(unpaired_output_file_name, "w")


#Using a series of dictionaries to sort reads
dict1 = dict()
dict2 = dict()
file_contents = open(input_fastq, 'r').read().split('\n')

for index in range(0, len(file_contents) -3, 4):
    header = file_contents[index].split('/')[0]
    seq = file_contents[index + 1]
    plus = file_contents[index + 2]
    qual = file_contents[index + 3]

    if header in dict1:
        dict2[header] = (seq, plus, qual)
        continue
    if not header in dict1:
        dict1[header] = (seq, plus, qual)


#Writing files based on presence in only one dictionary (unpaired) or both (paired)
for item in dict1:
    if not item in dict2:
        unpaired_output_file.write(item + '\n' + dict1[item][0] + '\n' + dict1[item][1] + '\n' + dict1[item][2] + '\n')
    else:
        paired_output_file.write(item + '/1' + '\n' + dict1[item][0] + '\n' + dict1[item][1] + '\n' + dict1[item][2] + '\n')
        paired_output_file.write(item + '/2' + '\n' + dict2[item][0] + '\n' + dict2[item][1] + '\n' + dict2[item][2] + '\n')


print('Output files generated for ' + OUT_PREFIX)
unpaired_output_file.close()
paired_output_file.close()

write_to_log(log_file, "process_read_pairs.py: Finished.")
