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
NOTE: to speed up the execution instead of 'SeqIO.parse()
low level fasta parser 'SimpleFastaParser' is used 
"""

import datetime
import os
import pandas as pd
from time import time
from Bio.SeqIO.FastaIO import SimpleFastaParser  # low level fast fasta parser
from ipywidgets import IntProgress
from IPython.display import display

NUCLEOTIDES = ["A", "T", "G", "C"]


def _get_coverage(input_file):
    """counts number of reads in fasta
    """
    read_counter = 0
    with open("./input_data/" + input_file) as in_handle:
        for record in SimpleFastaParser(in_handle):
            read_counter += 1
    # 1- without reference
    
    return read_counter - 1


def _get_ref_seq(input_file):
    """extracts first record from fasta
    which must be a reference
    """
    # low level parser returns tuple of id and sequence
    with open("./input_data/" + input_file) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            ref_seq = seq
            ref_seq_id = title
            break
        
    return ref_seq, ref_seq_id


def _split_gaps_from_pairs(ref_seq, seq):
    """cuts out pairs ("-", "-") from zipped
    reference and read sequences
    """
    pairs_ref_and_seq = list(zip(ref_seq, seq))
    pairs_gaps_cut = []
    for i in pairs_ref_and_seq:
        # skip gaps both in ref_seq and seq
        # as they mean insertion donwstream the alignment
        if i[0] == "-" and i[1] == "-":
            continue
        else:
            pairs_gaps_cut.append(i)

    return pairs_gaps_cut


def _split_ins_from_pairs(pairs_gaps_cut):
    """cuts out pairs ("-", some_nuc) from 
    zipped pairs"""
    pairs_ins_cut = []
    for i in pairs_gaps_cut:
        # skip insertions, as they shift true indexing
        if i[0] == "-" and i[1] in NUCLEOTIDES:
            continue
        else:
            pairs_ins_cut.append(i)

    return pairs_ins_cut


def _get_ref_array(ref_seq):
    """creates nested list of length equal 
    to reference sequence with
    all gaps cut off from reference
    """
    splitted_ref_seq = ref_seq.split("-")
    ref_true = "".join(splitted_ref_seq)
    ref_array = [[] for i in ref_true]

    return ref_array


def _get_indel_start_positions(seq):
    """creates list of indices which
    are indels start positions"""
    indel_positions = []
    indel_stop_index = 0

    for nuc_index in range(len(seq)):
        # as we can't change index of the next item in iteration
        # we skip nuc_index, until the indel stop index
        if nuc_index < indel_stop_index:
            continue

        # indel starts : "-"
        elif seq[nuc_index] == "-":
            indel_positions.append(nuc_index)
            # finding end of the indel
            for nuc_indel_index in range(nuc_index + 1, len(seq)):
                # it count like "A------" seq
                # considering all the gap till the end as true indel
                if seq[nuc_indel_index] == "-" and \
                        nuc_indel_index == len(seq) - 1:
                    indel_stop_index = nuc_indel_index + 1
                    break
                # indel continues: "-", skip
                elif seq[nuc_indel_index] == "-":
                    continue
                # we encounter nucleotide, so, indel ends
                # break this inner loop, and again find next indel
                elif seq[nuc_indel_index] in NUCLEOTIDES:
                    indel_stop_index = nuc_indel_index
                    break

    return indel_positions


def _count_indel_lens(seq_to_count_indels, indels_start_positions):
    """creates dictionary with positions as keys, and indels lengths 
    as values
    """
    pos_and_len = {k: None for k in indels_start_positions}
    indel_len_counter = 0

    for pos in indels_start_positions:

        for nuc_index in range(pos, len(seq_to_count_indels)):

            # catching single indel under the last nuc
            if seq_to_count_indels[nuc_index] == "-" and \
                    nuc_index == len(seq_to_count_indels) - 1:
                indel_len_counter += 1
                pos_and_len[pos] = indel_len_counter
                break

            elif seq_to_count_indels[nuc_index] == "-":
                indel_len_counter += 1

            # indel ends
            elif seq_to_count_indels[nuc_index] in NUCLEOTIDES:
                pos_and_len[pos] = indel_len_counter
                indel_len_counter = 0
                break

    return pos_and_len


def _seq_has_gaps_only(seq):

    if all(nuc == "-" for nuc in seq):
        result = True
    else:
        result = False

    return result


def _write_read_into_fasta(input_file, record, writ_mode):

    if not os.path.exists("./output_ids"):
        os.mkdir("output_ids")
    seq_file_id_name = "output_ids/" + "seq_ids_" + input_file

    with open(seq_file_id_name, writ_mode) as handle:
        handle.write(">" + record[0] + "\n" + record[1] + "\n")


def _count_dels(input_file, ref_seq, ref_seq_id, cov):

    total_dels = _get_ref_array(ref_seq)

    with open("./input_data/" + input_file) as in_handle:

        reads = SimpleFastaParser(in_handle)
        for record in reads:
            # skip first record as it's a reference
            if record[0] == ref_seq_id:
                # write reference into file with reads containing indels
                _write_read_into_fasta(input_file, record, "a")
                continue

            read = record[1]
            # pairs of nuc from ref and read without gaps ("-", "-")
            pairs_gaps_cut = _split_gaps_from_pairs(ref_seq, read)

            # this step is needed for dels
            # as insertions shift indexing forward
            pairs_ins_cut = _split_ins_from_pairs(pairs_gaps_cut)

            # take read, not ref - different from insertions counting
            seq_to_count = [x[1] for x in pairs_ins_cut]
            seq_to_count_dels = "".join(seq_to_count)

            # skip if read has noly gaps. it's not valid to count
            if _seq_has_gaps_only(seq_to_count_dels):
                continue

            dels_start_positions = _get_indel_start_positions(
                seq_to_count_dels)

            dels_pos_and_lens = _count_indel_lens(seq_to_count_dels,
                                                  dels_start_positions)

            # writing reads containing indels into a file
            if not all(value == None for key, value in dels_pos_and_lens.items()):
                _write_read_into_fasta(input_file, record, "a+")

            for key, value in dels_pos_and_lens.items():
                total_dels[key].append(value)

    return total_dels


def _count_ins(input_file, ref_seq, ref_seq_id, cov):

    total_ins = _get_ref_array(ref_seq)

    with open("./input_data/" + input_file) as in_handle:
        reads = SimpleFastaParser(in_handle)
        for record in reads:
            if record[0] == ref_seq_id:
                continue

            read = record[1]
            pairs_gaps_cut = _split_gaps_from_pairs(ref_seq, read)

            # take reference, not read. different from deletions counting
            seq_to_count = [x[0] for x in pairs_gaps_cut]
            seq_to_count_ins = "".join(seq_to_count)

            # skip if read has noly gaps. it's not valid to count
            if _seq_has_gaps_only(seq_to_count_ins):
                continue

            ins_start_positions = _get_indel_start_positions(seq_to_count_ins)

            ins_pos_and_lens = _count_indel_lens(seq_to_count_ins,
                                                 ins_start_positions)

            if not all(value == None for key, value in ins_pos_and_lens.items()):
                _write_read_into_fasta(input_file, record, "a+")

            correct_ins_pos_and_lens = {}
            len_correction = 0
            for key, value in ins_pos_and_lens.items():
                correct_ins_pos_and_lens[key - len_correction] = value
                len_correction += value

            for key, value in correct_ins_pos_and_lens.items():
                total_ins[key].append(value)

    return total_ins


def _create_df_cov(cov):
    df_cov = pd.DataFrame.from_dict({"coverage": cov}, orient='index').T
    return df_cov


def _create_df(ref_seq, total_indels):

    ref_no_gaps = ref_seq.split("-")
    ref = "".join(ref_no_gaps)
    columns = list(zip(ref, range(1, len(ref) + 1)))
    indel_df = pd.DataFrame(total_indels).T
    indel_df.columns = columns

    return indel_df


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
    
    ############################################
    # quick and dirty fix1 for reads that do not 
    # completely overlap the reference and s
    # easier to delete the data from the first and the last columns
    df_dels[df_dels.columns[0]] = None
    df_dels[df_dels.columns[-1]] = None
    df_ins[df_ins.columns[0]] = None
    #df_ins[df_ins.columns[-1]] = None # insertion under the last nuc may exist
    # end of fix1
    ############################################
    # fix2 for deletions, that start from the
    # 2nd nuc and later
    row_indexes = df_dels.index.values.tolist()
    column_indexes = df_dels.columns.values.tolist()
    ref_seq_len = len(column_indexes)
        
    for row_index in row_indexes:
        for col_index in range(0, len(column_indexes)):
            
            del_len = df_dels.iloc[row_index, col_index] # iloc as loc throws AssertionError somehow, have no time to deep into details
            if del_len == (ref_seq_len - col_index):
                df_dels.iloc[row_index, col_index] = None
    # end of fix2
    ###############################################
    # fix3 : clean the df after, deleting empty rows
    df_dels.dropna(axis=0, how="all", inplace=True)
    df_ins.dropna(axis=0, how="all", inplace=True)
    # end of fix3
    #############################################
    
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
        input_files = [f for f in input_files if f.rsplit(
            ".", 1)[-1] in extensions]

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
                      Now this file is skipped: {1}
                       """.format(f, permerr))
            except IndexError as inderr:
                print("something wrong with {0}, error: {1}".format(f, inderr))
            try:
                df_dels = _create_df(ref_seq, total_dels)
                df_ins = _create_df(ref_seq, total_ins)
                df_cov = _create_df_cov(cov)

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
