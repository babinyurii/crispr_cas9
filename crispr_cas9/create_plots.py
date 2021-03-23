# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:21:42 2019

@author: babin
"""

import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import datetime
from time import time
from ipywidgets import IntProgress
from IPython.display import display
sns.set(font_scale=1.5)


def _get_current_time():
    """returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')
    return time_stamp


def _create_bars(file_name, indel_matrix, indel):
    """percent reads count 
    based on coverage
    """
    raw_del_count = indel_matrix.sum(axis=0)
    fig_del_percent = plt.figure(figsize=(20, 7))
    #perc_del_count = raw_del_count / cov * 100

    plt.bar(raw_del_count.index, raw_del_count)

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
    #largest_indel_num = indel_matrix.max()
    # turning into int class numpy float, otherwise deprecation warning is raised
    #largest_indel_num = int(largest_indel_num.max())

    fig = plt.figure(figsize=(25, 12))
    ax = sns.heatmap(indel_matrix, cmap=cmap_short, annot=False,
                     # linewidths=1, linecolor='grey',
                     cbar_kws={'label': indel + " percent"})

    # setting font size of the colormap label
    ax.figure.axes[-1].yaxis.label.set_size(25)

    plt.yticks(list(range(1, len(indel_matrix.index) + 1)),
               list(reversed(range(1, len(indel_matrix.index) + 1))),
               fontsize=12, horizontalalignment="left")
    plt.xlabel("nucleotide position", fontsize=20)
    plt.ylabel("length of indel", fontsize=20)
    plt.title("percent and length of " + indel, fontsize=25)
    # plt.show() # comment this line to save figure, otherwise it'll save blank file
    fig.savefig("./output_plots/" + file_name.rsplit(".", 1)
                [0] + "_heatmap_" + indel + ".png")
    plt.close(fig)  # comment the line to show the figure in the jupyter


def _show_report(total_time, file_counter):
    """prints out very brief report 
    """
    hours = total_time // 3600
    minutes = (total_time % 3600) // 60
    seconds = total_time % 60

    print("""
    
    file processed: {0}
    time taken: {1} hours {2} minutes {3} seconds
    
    the results are in the folder 'output_plots'
    
    ---------------
    ... job finished at {4}
    ---------------
    
    """.format(file_counter, hours, minutes, int(seconds), _get_current_time()))


def main():
    """main functions
    """
    start_time = time()
    file_counter = 0

    print("""
    ---------------
    job started at {0} ...
    ---------------
    """.format(_get_current_time()))

    if not os.path.exists("./output_plots"):
        os.mkdir("output_plots")

    if os.path.exists("./output_matrices"):
        input_sheets = os.listdir("./output_matrices")
        input_sheets = [
            f for f in input_sheets if f.rsplit(".", 1)[-1] == "xlsx"]

        num_files = len(input_sheets)
        progress_bar = IntProgress(min=0, max=num_files, bar_style='success')
        display(progress_bar)

        for f in input_sheets:
            try:
                file_name = f.rsplit(".", 1)[0]
                file_content = file_name.split("_")[-2]
                indel_type = file_name.split("_")[-1]

                indel_matrix = pd.read_excel("./output_matrices/" + f)
                indel_matrix.index = indel_matrix["indel length"]
                indel_matrix = indel_matrix.drop(columns=["indel length"])

                if file_content == "count":
                    _create_bars(f, indel_matrix, indel_type)
                elif file_content == "percent":
                    _create_heatmap(f, indel_matrix, indel_type)

            except PermissionError as perr:
                print("""
                              warning: perhaps the excel spreadsheet is open in excel. 
                              please, close it and rerun the scrips. File {0} will be skipped.
                              error: {1}
                              """.format(f, perr))

            progress_bar.value += 1
            file_counter += 1
        # finally after all the files were iterated through
        else:
            finish_time = time()
            total_time = finish_time - start_time
            _show_report(total_time, file_counter)

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
    main()
