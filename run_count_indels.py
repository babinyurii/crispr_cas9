# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 17:02:21 2018
@author: babin
"""
""" run_count_indels.py
this script count exact number of indels in "fasta" file
and writes results into the excel spreadsheet.
it takes 'fasta' file imported from Geneious as an input
'fasta' file must be imported from Geneious with the following parameters:
    - wrap sequence lines every 80 chars: No
    - replace spaces in sequence names with undescores: Yes
    - include sequence description: Yes
    - export sequence in: upper case
    - export missing ends of alignments as: Ns
to use the script:
    1. create folder 'input_data' in the current directory and paste 'fastas' in it
    2. import script into the current directory using Jupyter or whatever and run it:
        from run_count_indels import run_count_indels
        run_count_indels()
    or put the script into the current directory and run via shell:
        python run_count_indels.py
"""


import datetime
import os
import pandas as pd
from Bio import SeqIO
from collections import OrderedDict
from time import time



FILE_COUNTER = 0

def _count_indels(input_file):
    """counts deletions and insertions,
    collects them into nested lists.
    each nested list in the array
    corresponds to the position 
    in the mapping ('raw', with gaps
    at the insertion sites)
    """
    nucleotides = ["A", "T", "G", "C"]

    for seq in SeqIO.parse("./input_data/" + input_file, "fasta"):
        ref_seq = seq
        break

    total_deletions = [[] for x in range(len(ref_seq))]
    total_insertions = [[] for x in range(len(ref_seq))]

    reads = SeqIO.parse("./input_data/" + input_file, "fasta")

    for read in reads:
        if read.id == ref_seq.id:
            continue

        start_nuc = 0

        for nuc in range(len(read)):
            if nuc < start_nuc:
                continue

            elif ref_seq[nuc] in nucleotides and read[nuc] in nucleotides:
                continue

            elif ref_seq[nuc] == read[nuc]:
                continue

            # counting deletion length
            elif ref_seq[nuc] in nucleotides and read[nuc] == "-":
                deletions_counter = 0

                for nuc_del in range(nuc + 1, len(read)):
                    if read[nuc_del] == "-" and ref_seq[nuc_del] in nucleotides:
                        deletions_counter += 1
                    elif read[nuc_del] == "-" and ref_seq[nuc_del] == "-":
                        continue
                    elif read[nuc_del] in nucleotides or read[nuc_del] in nucleotides and ref_seq[nuc_del] in nucleotides:
                        start_nuc = nuc_del
                        total_deletions[nuc].append(deletions_counter + 1)
                        break

            # counting insertion length
            elif ref_seq[nuc] == "-" and read[nuc] in nucleotides:
                insertion_counter = 0

                for nuc_ins in range(nuc + 1, len(read)):
                    if ref_seq[nuc_ins] == "-" and read[nuc_ins] in nucleotides:
                        insertion_counter += 1
                    elif read[nuc_ins] == "-" and ref_seq[nuc_ins] == "-":
                        continue
                    elif ref_seq[nuc_ins] in nucleotides or read[nuc_del] in nucleotides and ref_seq[nuc_ins] in nucleotides:
                        start_nuc = nuc_ins
                        total_insertions[nuc].append(insertion_counter + 1)
                        break

    file_name = input_file.rsplit(".", 1)[0]
    _create_df(ref_seq, total_deletions, total_insertions, file_name)


def _create_df(ref_seq, total_deletions, total_insertions, file_name):
    """processes deletion and insertions lists.
    corrects nucleotides indices 
    by original reference sequence
    (indexing without gaps)
    creates pandas DataFrames, 
    writes them into excel spreadsheets
    """
    ref_raw = str(ref_seq.seq)  # converting SeqRecord obj into python str type
    # adding indices to the nucleotides, by ref with gaps
    ref_indexed = list(zip(ref_raw, range(len(ref_raw))))
    ref_dels = [x for x in ref_indexed if x[0] != "-"]  # removing gaps

    # creating 'real' string which were due to mapping
    ref = ref_raw.split("-")
    ref = "".join(ref)

    # columns for data frames, indices start with 1
    columns = list(zip(ref, range(1, len(ref) + 1)))

    # collecting deletions
    # index of nucleotides and array with gaps correspond to the same position
    container_del = OrderedDict()
    for item in ref_dels:
        container_del[item] = total_deletions[item[1]]

    # collecting insertions
    container_ins = OrderedDict()
    for item in ref_indexed:
        container_ins[item] = total_insertions[item[1]]
    # taking only the positions with insertions
    container_ins_correct = OrderedDict()
    insertion_counter = 0

    for key, value in container_ins.items():
        # if so, the position (in the ref_seq with gaps) has insertion under it
        if key[0] == "-" and value != []:
            # adding 1 as due to the insertions and gaps real indices are shifted forward
            insertion_counter += 1
            # add insertion under the real index (in ref_seq with no gaps)
            container_ins_correct[key[1] - insertion_counter] = value
        elif key[0] == "-" and value == []:
            insertion_counter += 1  # add 1 to counter, as gaps shift it as well

    # adding lost indices from the ref with no gaps
    for ind in range(len(ref)):
        if ind not in container_ins_correct.keys():  # if index not in keys
            # add it to the array with empty list
            container_ins_correct[ind] = []
    # sort by added indices
    container_ins_correct = OrderedDict(
        sorted(container_ins_correct.items(), key=lambda x: x[0]))

    df_dels = pd.DataFrame.from_dict(container_del, orient='index').T
    df_ins = pd.DataFrame.from_dict(container_ins_correct, orient='index').T
    df_dels.columns = columns
    df_ins.columns = columns

    if not os.path.exists("./output_indels"):
        os.mkdir("output_indels")

    writer = pd.ExcelWriter("./output_indels/" + file_name + '.xlsx')
    df_dels.to_excel(writer, "deletions")
    df_ins.to_excel(writer, "insertions")
    
    try:
        writer.save()
        global FILE_COUNTER
        FILE_COUNTER += 1
        _get_current_time()
        print("file '{0}' processing finished at: {1} \n".format(
        file_name, _get_current_time()))
    except PermissionError as e:
        print(
            """
        Houston, we have a problem...
        ---------
        It seems that excel spreadsheet '{0}' that was going to be rewritten 
        is open in excel. To process the corresponding 'fasta' file, rerun the script
        Now this file is skipped
        ---------
        """.format(file_name))


def _get_current_time():
    """just returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')
    return time_stamp


def _show_report(total_time):
    """prints out very brief report 
    """
    global FILE_COUNTER
    hours = total_time // 3600
    minutes = (total_time % 3600) // 60
    seconds = total_time % 60
    
    total_files = FILE_COUNTER

    print("""
    
    file processed: {0}
    time taken: {1} hours {2} minutes {3} seconds
    
    the results are in the folder 'output_indels'
    
    ---------------
    ... job finished at {4}
    ---------------
    
    """.format(total_files, hours, minutes, int(seconds), _get_current_time()))
    FILE_COUNTER = 0
    
    

def run_count_indels():
    """main functions
    """
    start_time = time()
    
    print("""
    ---------------
    job started at {0} ...
    ---------------
    """.format(_get_current_time()))
    
    extensions = ["fasta", "fa", "fas"]
    if os.path.exists("./input_data"):
        input_files = os.listdir("./input_data")

        for f in input_files:
            if os.path.isdir("./input_data/" + f):
                print("warning: please, don't put folders into the 'input_data' folder.\
                I believe this is a folder: '{0}':".format(f) + "\n")
            elif f.rsplit(".", 1)[-1] in extensions:
                print("processing file '{}'".format(f))
                _count_indels(f)
            else:
                print(
                    "warning: item '{0}' isn't a fasta file. it won't be processed\n".format(f))
        else:
            finish_time = time()
            total_time = finish_time - start_time
            _show_report(total_time)
    else:
        os.mkdir("input_data")
        print(
            """
        Houston, we have a problem...
        --------
        folder 'input_data' doesn't exist in the current directory 
        or may be you've created it but misspelled its name.
        Anyway, it has just beencreated by this script.
        Paste your 'fasta' files into the folder 'input_data' and run this script again.
        --------
        """
        )


if __name__ == "__main__":
    run_count_indels()