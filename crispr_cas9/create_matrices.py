# -*- coding: utf-8 -*-

import datetime
import numpy as np
import pandas as pd
import os
from collections import OrderedDict
from time import time
from xlrd import XLRDError
from ipywidgets import IntProgress
from IPython.display import display

def _df_test(file_name):
    df_test = pd.read_excel("./output_indels/" + file_name)
    print(df_test)

def _get_coverage_from_excel(file_name):
    """takes out coverage value from
    excel file.
    """
    df_cov = pd.read_excel("./output_indels/" + file_name, "coverage")
    cov = int(df_cov["coverage"])
    return cov


def _get_current_time():
    """returns time stamp
    """
    time_stamp = datetime.datetime.fromtimestamp(
        time()).strftime('%Y-%m-%d %H:%M:%S')
    return time_stamp


def _save_indel_count_matrix(file_name, indel_matrix, indel):

    writer = pd.ExcelWriter("./output_matrices/" + file_name.rsplit(".", 1)[0]
                            + '_matrix_count_' + indel + '.xlsx')
    indel_matrix.index.name = "indel length"
    
    
    indel_matrix.drop(indel_matrix.columns[indel_matrix.columns.str.contains("Unnamed", case=False)],
                     axis=1, inplace=True)
    
    indel_matrix.to_excel(writer)
    writer.save() #FIXED

def _save_indel_percent_matrix(file_name, indel_matrix_percent, indel):

    writer = pd.ExcelWriter("./output_matrices/" + file_name.rsplit(".", 1)[0]
                            + '_matrix_percent_' + indel + '.xlsx')
    indel_matrix_percent.index.name = "indel length"
    indel_matrix_percent.drop(indel_matrix_percent.columns[indel_matrix_percent.columns.str.contains("Unnamed", case=False)],
                     axis=1, inplace=True)
    indel_matrix_percent.to_excel(writer)
    writer.save() #FIXED

def _create_matrix(file_name, indel):
    """
    creates matrix (length of indel, position) = number of such indels
    """
    # creating matrix of rows:length of indel, columns:position in the sequence
    df_indels = pd.read_excel("./output_indels/" + file_name, sheet_name=indel, index_col=0)

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
    max_col_ind = int(max(indel_matrix.columns))
    for i in range(1, max_col_ind):
        if i not in indel_matrix.columns:
            # adding missing value ('indel length') to column with empty values
            indel_matrix[i] = np.nan

    indel_matrix = indel_matrix.T
    indel_matrix.sort_index(axis=0, ascending=False, inplace=True)
    indel_matrix.columns = df_indels.columns
    indel_matrix = indel_matrix.fillna(value=0)

    return indel_matrix


def _create_matrix_percent(file_name, indel_matrix, cov, indel):
    """converts raw count value in matrix 
    into percentage based on coverage value"""

    for ind in indel_matrix.index:
        for col in indel_matrix.columns:
            indel_matrix.loc[ind, col] = indel_matrix.loc[ind, col] / cov * 100

    return indel_matrix


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
    """main functions
    """
    start_time = time()
    file_counter = 0

    print("""
    ---------------
    job started at {0} ...
    ---------------
    """.format(_get_current_time()))

    indel_type = ["deletions", "insertions"]

    if not os.path.exists("./output_matrices"):
        os.mkdir("output_matrices")

    if os.path.exists("./output_indels"):
        input_sheets = os.listdir("./output_indels")
        input_sheets = [
            f for f in input_sheets if f.rsplit(".", 1)[-1] == "xlsx"]

        num_files = len(input_sheets)
        progress_bar = IntProgress(min=0, max=num_files, bar_style='success')
        display(progress_bar)
        for f in input_sheets:
            try:
                cov = _get_coverage_from_excel(f)
            except (PermissionError, KeyError, XLRDError) as e:
                print("warning: coverage can't be derived from the file {0}, error {1}".format(f, e))
                continue

            for indel in indel_type:
                try:
                    indel_matrix = _create_matrix(f, indel)
                    _save_indel_count_matrix(f, indel_matrix, indel)

                    matrix_percent = _create_matrix_percent(f, indel_matrix,
                                                            cov, indel)
                    _save_indel_percent_matrix(f, matrix_percent, indel)

                except PermissionError as perr:
                    print("""
                          warning: perhaps the excel spreadsheet is open in excel. 
                          please, close it and rerun the scrips. File {0} will be skipped.
                          error: {1}
                          """.format(f, perr))
                except KeyError as kerr:
                    print("""
                          warning: I believe the file {0} somehow is not valid. 
                          It wil be skipped. 
                          It may be an empty excel spreadsheet, or broken file with results. 
                          Please, check out, whether is has three tabs:
                          'deletions', 'insertions', 'coverage'
                          error: {1}
                          """.format(f, kerr))
                except XLRDError as xlrerr:
                    print("warning: no tab. error: {0}".format(xlrerr))
                except ValueError as valerr:
                    print("""
                          warning: perhaps the excel tab "{3}" has no data.
                          the {0} for file {1} won't be written.
                          error: {2}.
                          """.format(indel, f, valerr, indel))
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
