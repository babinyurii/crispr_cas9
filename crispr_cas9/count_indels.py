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
IMPORTANT NOTES:
    - first record in fasta file must be your reference sequence
    - reference sequence must start and end with nucleotides, not gaps
    f.e.: 'A-A---A---A' is the right string, 
    but: '----A---AAAA---' is not valid
    if reference sequence starts or ends with gap, the scrips will throw: 
    'ValueError: Length mismatch: Expected axis has N elements, 
    new values have n elements'. This is the pandas exception concerning indices.

to use the script:
    1. create folder 'input_data' in the current directory and paste 'fasta' files into it
    2. import script into the current directory using Jupyter or whatever and run it:
        1) from crispr_count_indels import run_count_indels
        run_count_indels()
        2) %run crispr_count_indels.py
    or put the script into the current directory and run via shell:
    python run_count_indels.py

NOTE: to speed up the exection instead of 'SeqIO.parse()
low level fasta parser 'SimpleFastaParser' is used 
"""


import datetime
import os
import pandas as pd
from collections import OrderedDict
from time import time
from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser
from ipywidgets import IntProgress
from IPython.display import display




def _get_coverage(input_file):
    read_counter = 0
    with open("./input_data/" + input_file) as in_handle:
        for record in SimpleFastaParser(in_handle):
            read_counter += 1
    # 1- without reference        
    return read_counter - 1

def _get_ref_seq(input_file):
    # low level parser returns tuple of id and sequence
    with open("./input_data/" + input_file) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            ref_seq = seq
            ref_seq_id = title
            break
    return ref_seq, ref_seq_id


    
# crazy _count_indels() is the result of fast legacy rewriting
# TODO
# 4. rename vars. f.e. nuc is not nucleotide, its nucleotide position counter actually
    
def _count_dels(input_file, ref_seq, ref_seq_id, cov):
    """counts deletions and insertions, collects them into nested lists.
    each nested list in the array corresponds to the position 
    in the mapping ('raw', with gaps at the insertion sites)
    """
    nucleotides = ["A", "T", "G", "C"]
    total_deletions = [[] for x in range(len(ref_seq))]

    with open("./input_data/" + input_file) as in_handle:
        reads = SimpleFastaParser(in_handle)

        for record in reads:
            read = record[1]
            start_nuc = 0
            for nuc in range(len(read)):
                # if the nucleotide to slice from is the last
                # and the indel is over, we must stop
                # as the stopping condition was accomplished:
                if start_nuc == len(ref_seq) - 1:
                    break

                elif nuc < start_nuc:
                    continue
                # handling single deletion under the last nucleotide
                elif ref_seq[nuc] in nucleotides and read[nuc] == "-" and \
                    nuc == len(read) - 1:
                        total_deletions[nuc].append(1)
                        break
                # counting deletion length
                elif ref_seq[nuc] in nucleotides and read[nuc] == "-":
                    deletions_counter = 1

                    for nuc_del in range(nuc + 1, len(read)):  
                        # deletions end var 0: the final nucleotide is deleted
                        # check up length for verify
                        if ref_seq[nuc_del] in nucleotides and read[nuc_del] == "-" and \
                        nuc_del == len(ref_seq) - 1:
                            start_nuc = nuc_del
                            total_deletions[nuc].append(deletions_counter + 1)
                            break
                        # deletion goes on
                        elif ref_seq[nuc_del] in nucleotides and read[nuc_del] == "-":
                            deletions_counter += 1
                        # skip if gap is both in ref and read
                        # as it means insertion somewhere downstream the file
                        elif ref_seq[nuc_del] == "-" and read[nuc_del] == "-":
                            continue
                        # deletion ends var 1 : both ref and read are in nucleotides
                        elif ref_seq[nuc_del] in nucleotides and read[nuc_del] in nucleotides:
                            start_nuc = nuc_del 
                            total_deletions[nuc].append(deletions_counter)
                            break
                        # deletion ends var 2 : ref seq has a gap and read in nucleotides,
                        # which is an insertion, not a deletion
                        elif ref_seq[nuc_del] == "-" and read[nuc_del] in nucleotides:
                            start_nuc = nuc_del 
                            total_deletions[nuc].append(deletions_counter)
                            break
    return total_deletions


def _count_ins(input_file, ref_seq, ref_seq_id, cov):
    
    nucleotides = ["A", "T", "G", "C"]
    total_insertions = [[] for x in range(len(ref_seq))]

    with open("./input_data/" + input_file) as in_handle:
        reads = SimpleFastaParser(in_handle)

        for record in reads:
            read = record[1]
            start_nuc = 0

            for nuc in range(len(read)):
                
                # if the nucleotide to slice from is the last
                # and the indel is over, we must stop
                if start_nuc == len(ref_seq) - 1:
                    break
                
                elif nuc < start_nuc:
                    continue
                
                # counting insertion length
                elif ref_seq[nuc] == "-" and read[nuc] in nucleotides:
                    insertion_counter = 0

                    for nuc_ins in range(nuc + 1, len(read)):
                        # insertion goes on
                        if ref_seq[nuc_ins] == "-" and read[nuc_ins] in nucleotides:
                            insertion_counter += 1
                        # skip if gap is both in ref and read
                        #it means insertions somewhere downstream the file
                        elif read[nuc_ins] == "-" and ref_seq[nuc_ins] == "-":
                            continue
                        # insertions end var 1 : both ref and read are in nucleotides
                        elif ref_seq[nuc_ins] in nucleotides and read[nuc_ins] in nucleotides:
                            start_nuc = nuc_ins
                            total_insertions[nuc].append(insertion_counter + 1)
                            break
                        # insertion end var 2 : ref seq is in nucleotides, 
                        # and read has a gap, which is a deletion, not a insertion
                        elif ref_seq[nuc_ins] in nucleotides and read[nuc_ins] == "-":
                            start_nuc = nuc_ins
                            total_insertions[nuc].append(insertion_counter + 1)
                            break

    return total_insertions


def _create_df_cov(cov):
    df_cov = pd.DataFrame.from_dict({"coverage": cov}, orient='index').T
    return df_cov

    
def _create_df_dels(ref_seq, total_deletions, file_name, cov):
    """processes deletion and insertions lists.
    corrects nucleotides indices 
    by original reference sequence
    (indexing without gaps)
    creates pandas DataFrames, 
    writes them into excel spreadsheets
    """
    ref_raw = ref_seq
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
    
    df_dels = pd.DataFrame.from_dict(container_del, orient='index').T
    df_dels.columns = columns
    
    return df_dels

# logic here is rather awkward, as insertion position should be correced 
# because of the shift
def _create_df_ins(ref_seq, total_insertions, file_name, cov):
    """processes deletion and insertions lists.
    corrects nucleotides indices 
    by original reference sequence
    (indexing without gaps)
    creates pandas DataFrames, 
    writes them into excel spreadsheets
    """
    ref_raw = ref_seq
    # adding indices to the nucleotides, by ref with gaps
    ref_indexed = list(zip(ref_raw, range(len(ref_raw))))

    # creating 'real' string which were due to mapping
    ref = ref_raw.split("-")
    ref = "".join(ref)

    # columns for data frames, indices start with 1
    columns = list(zip(ref, range(1, len(ref) + 1)))

    container_ins = OrderedDict()
    for item in ref_indexed:
        container_ins[item] = total_insertions[item[1]]
    # taking only the positions with insertions
    container_ins_correct = OrderedDict()
    insertion_counter = 0

    for key, value in container_ins.items():
        # if so, the position (in the ref_seq with gaps) has insertion under it
        if key[0] == "-" and value != []:
            
            # add insertion under the real index (in ref_seq with no gaps)
            container_ins_correct[key[1] - insertion_counter] = value
            # adding 1 as due to the insertions which 'move' real indices forward
            # incrementing after as first insertion is under the 'real', without 
            # gaps index
            insertion_counter += 1
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
    
    df_ins = pd.DataFrame.from_dict(container_ins_correct, orient='index').T
    df_ins.columns = columns
    
    return df_ins



def _get_current_time():
    """just returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')
    return time_stamp


def _save_df(df_dels, df_ins, df_cov, file_name):
    """saves dataframes into an excel spreadsheet
    """
    if not os.path.exists("./output_indels"):
        os.mkdir("output_indels")

    file_name = file_name.rsplit(".", 1)[0]
    writer = pd.ExcelWriter("./output_indels/" + file_name + '.xlsx')
    df_dels.to_excel(writer, "deletions")
    df_ins.to_excel(writer, "insertions")
    df_cov.to_excel(writer, "coverage")

    writer.save()
   

def _show_report(total_time, file_counter):
    """prints out very brief report 
    """

    hours = total_time // 3600
    minutes = (total_time % 3600) // 60
    seconds = total_time % 60

    print("""
    
    file processed: {0}
    time taken: {1} hours {2} minutes {3} seconds
    
    the results are in the folder 'output_indels'
    
    ---------------
    ... job finished at {4}
    ---------------
    
    """.format(file_counter, hours, minutes, int(seconds), _get_current_time()))


def main():
    """main function
    """
    start_time = time()
    file_counter = 0
    
    print("""
    ---------------
    job started at {0} ...
    ---------------
    """.format(_get_current_time()))
    
    extensions = ["fasta", "fa", "fas"]
    if os.path.exists("./input_data"):
        input_files = os.listdir("./input_data")
        input_files = [f for f in input_files if f.rsplit(".", 1)[-1] in extensions]
        
        num_files = len(input_files)
        progress_bar = IntProgress(min=0, max=num_files, bar_style='success')
        display(progress_bar)
        
        for f in input_files:            
            try:
                ref_seq, ref_seq_id = _get_ref_seq(f)
                cov = _get_coverage(f)
                total_dels = _count_dels(f, ref_seq, ref_seq_id, cov)
                total_ins = _count_ins(f, ref_seq, ref_seq_id, cov)
            except PermissionError as permerr:
                print("""
                      warning: It seems that the file'{0}' is open in some other 
                      application which doesn't allow it to be processed
                      Now this file is skipped
                       """.format(f))
            try:
                df_dels = _create_df_dels(ref_seq, total_dels, f, cov)
                df_ins = _create_df_ins(ref_seq, total_ins, f, cov)
                df_cov = _create_df_cov(cov)
                
                
                #df_dels, df_ins, df_cov = _create_df(ref_seq, total_dels,
                                                # total_ins, f, cov)
            except ValueError as valerr:
                print("""
                      warning: it seems that the file {0}
                      has no reference sequence as the first record. 
                      please, check it out
                      the file will be skipped. error: {1}
                      """.format(f, valerr))
            try:
                _save_df(df_dels, df_ins, df_cov, f)
            except PermissionError as permerr:
                print("""
                      warning: It seems that excel spreadsheet '{0}' is open in excel.
                      To process the corresponding 'fasta' file, rerun the script
                      now this file is skipped. error: {1}
                      """.format(f, permerr))
            
            progress_bar.value += 1
            file_counter += 1 
                
        finish_time = time()
        total_time = finish_time - start_time
        _show_report(total_time, file_counter)
    
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
    main()
