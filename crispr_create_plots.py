# -*- coding: utf-8 -*-

""" crispr_create_plots.py
this script creates bar charts and heatmaps using 
output excel spreadsheets made by 'crispr_count_indels.py'.
bar charts show percent of reads with indels irrespective of the indel length
and give general overview of the data
heatmaps show exact number of indels and their length

to use the script:
    1. create folder 'input_data' in the current directory and paste 'fasta' files into it
    2. import script into the current directory using Jupyter run it:
        a) from run_count_indels import run_count_indels
        run_count_indels()
        b) %run crispr_create_plots.py
    or put the script into the current directory and run via shell:
    python run_create_plots.py
"""

import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns
sns.set()
from collections import OrderedDict
from time import time


FILE_COUNTER = 0


def _get_coverage(file_name):
    """takes out coverage value from
    excel file.
    """
    df_cov = pd.read_excel(file_name, 2)
    cov = int(df_cov["coverage"])
    return cov


def _get_current_time():
    """returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')
    return time_stamp


def _create_matrix(file_name):
    """
    creates matrix (length of indel, position) = number of such indels
    """

    cov = _get_coverage("./output_indels/" + file_name)

    indel_type = ["deletions", "insertions"]

    # creating matrix of rows:length of indel, columns:position in the sequence
    for indel in indel_type:
        df_indels = pd.read_excel(
            "./output_indels/" + file_name, sheet_name=indel)

        container = []

        for col in df_indels.columns:
            indels_total = df_indels.loc[:, col].tolist()
            lengths = list(set(indels_total))
            lengths = (x for x in lengths if x == x)  # deleting 'nan'
            d = OrderedDict()

            for length in lengths:
                d[length] = indels_total.count(length)
            container.append(d)

        indel_matrix = pd.DataFrame(container)

        # adding missing values to df columns with empty values
        # as we don't have all lengths of indels
        # it's for accurate plots
        try:
            max_col_ind = int(max(indel_matrix.columns))
            for i in range(1, max_col_ind):
                if i not in indel_matrix.columns:
                    # adding missing value ('indel length') to column with empty values
                    indel_matrix[i] = np.nan

            indel_matrix = indel_matrix.T
            indel_matrix.sort_index(axis=0, ascending=False, inplace=True)
            indel_matrix.columns = df_indels.columns

            _create_bars_percent(file_name, indel_matrix, cov, indel)
            _create_heatmap(file_name, indel_matrix, indel)

        # trying to catch empty excel tabs
        except ValueError as e:
            print("""warning: perhaps the excel tab is empty. 
            the {0} for file {1} won't be written.
            error: {2}""".format(indel, file_name, e))

        # return indel_matrix


def _create_bars_percent(file_name, indel_matrix, cov, indel):
    """percent reads count 
    based on coverage
    """
    raw_del_count = indel_matrix.sum(axis=0)
    fig_del_percent = plt.figure(figsize=(20, 7))
    perc_del_count = raw_del_count / cov * 100

    plt.bar(perc_del_count.index, perc_del_count)
    plt.xticks(list(range(len(raw_del_count))),
               list(range(1, len(raw_del_count)+1)))

    plt.title("percent of " + indel, fontsize=25)
    plt.xlabel("nucleotide position", fontsize=25)
    plt.ylabel("percent", fontsize=25)

    if not os.path.exists("./output_plots"):
        os.mkdir("output_plots")
    # plt.show() # comment to save figure, otherwise it'll save blank file
    fig_del_percent.savefig(
        "./output_plots/" + file_name.rsplit(".", 1)[0] + "_bars_" + indel + ".png")
    plt.close()  # comment the line to show the figure in the jupyter or wherever


def _create_heatmap(file_name, indel_matrix, indel):

    cmap = plt.cm.jet  # defining colormap
    cmaplist = [cmap(i) for i in range(cmap.N)]  # list of colors from jet map
    cmap_short = cmaplist[130:]
    cmap_short[0] = (.5, .5, .5, 1.0)  # first entry changed to be grey

    # here's the max value of indels in the current df
    # it's used to denote the largest value in colormap and split colormap into segments
    largest_indel_num = indel_matrix.max()
    # turning into int class numpy float, otherwise deprecation warning is raised
    largest_indel_num = int(largest_indel_num.max())

    fig = plt.figure(figsize=(25, 8))
    #ax = plt.axes()

    sns.heatmap(indel_matrix, cmap=cmap_short, annot=True,
                cbar_kws={'label': indel + " count"})
    plt.yticks(list(range(1, len(indel_matrix.index) + 1)),
               list(reversed(range(1, len(indel_matrix.index) + 1))),
               fontsize=12, horizontalalignment="left")

    plt.xlabel("nucleotide position", fontsize=20)
    plt.ylabel("length of indel", fontsize=20)

    plt.title("number and length of " + indel, fontsize=25)
    #plt.xticks(list(range(1, len(matrix.columns) + 1)), list(range(1, len(matrix.columns) + 1)), fontsize=15)
    # plt.show() # comment to save figure, otherwise it'll save blank file
    fig.savefig("./output_plots/" + file_name.rsplit(".", 1)
                [0] + "_heatmap_" + indel + ".png")
    plt.close()  # comment the line to show the figure in the jupyter or wherever


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
    
    the results are in the folder 'output_plots'
    
    ---------------
    ... job finished at {4}
    ---------------
    
    """.format(total_files, hours, minutes, int(seconds), _get_current_time()))
    FILE_COUNTER = 0


def run_create_plots():
    """main functions
    """
    start_time = time()

    print("""
    ---------------
    job started at {0} ...
    ---------------
    """.format(_get_current_time()))

    if os.path.exists("./output_indels"):
        input_sheets = os.listdir("./output_indels")

        for f in input_sheets:
            if os.path.isdir("./output_indels/" + f):
                print("""\nwarning: please, don't put folders into the 'output_indels' folder
                I believe this is a folder: '{0}'""".format(f) + "\n")
            elif f.rsplit(".", 1)[-1] == "xlsx":
                print("processing file '{}'".format(f))

                # if spreadsheet is open in excel skip it for now
                try:
                    _create_matrix(f)
                    global FILE_COUNTER
                    FILE_COUNTER += 1
                    print("file '{0}' processing finished at: {1} \n".format(
                        f, _get_current_time()))
                except PermissionError as perr:
                    print("""warning: perhaps the excel spreadsheet is open in excel. 
                    please, close it and rerun the scrips. File {0} will be skipped.
                    error: {1}""".format(f, perr))
                except KeyError as kerr:
                    print("""warning: I believe the file {0} somehow is not valid. It wil be skipped. 
                    It may be an empty excel spreadsheet, or broken file with results. 
                    Please, check out, whether is has three tabs: 'deletions', 'insertions', 'coverage'.
                    error: {1}""".format(f, kerr))
            else:
                print(
                    "\nwarning: item '{0}' isn't a fasta file. it won't be processed\n".format(f))

        # finally after all the files were iterated through
        else:
            finish_time = time()
            total_time = finish_time - start_time
            _show_report(total_time)
    else:
        print(
            """
        Houston, we have a problem...
        --------
        folder 'output_indels' doesn't exist in the current directory 
        or may be you've renamed it but misspelled its name.
        Please, check it out.
        --------
        """
        )


if __name__ == "__main__":
    run_create_plots()
