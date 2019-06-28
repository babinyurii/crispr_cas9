# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:28:01 2019

@author: babin
"""

import sys
sys.path.append("..")

from crispr_cas9 import count_indels
from crispr_cas9 import create_matrices


def test_create_plots_create_matrix_deletion_raw_count():
    """
    """
    #ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("15_dels_ladder_for_matrix.fasta")
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("15_dels_ladder_for_matrix.fasta")
    cov = count_indels._get_coverage("15_dels_ladder_for_matrix.fasta")
    total_dels = count_indels._count_dels("15_dels_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("15_dels_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    
    
    df_dels = count_indels._create_df(ref_seq, total_dels)
    df_ins = count_indels._create_df(ref_seq, total_ins)
    df_cov = count_indels._create_df_cov(cov)
    
    count_indels._save_df(df_dels, df_ins, df_cov, "15_dels_ladder_for_matrix.fasta")
    
    deletion_matrix = create_matrices._create_matrix("15_dels_ladder_for_matrix.xlsx", "deletions")
    try:
        insertion_matrix = create_matrices._create_matrix("15_dels_ladder_for_matrix.xlsx", "insertions")
    except ValueError:
        print("no insertion matrix")
        
    
    # non empty values, counted indels
    assert deletion_matrix.iloc[4, 0] == 2
    assert deletion_matrix.iloc[3, 1] == 2
    assert deletion_matrix.iloc[2, 2] == 2
    assert deletion_matrix.iloc[1, 3] == 2
    assert deletion_matrix.iloc[0, 4] == 2
    assert deletion_matrix.iloc[0, 5] == 2
    assert deletion_matrix.iloc[1, 6] == 1
    assert deletion_matrix.iloc[2, 7] == 1
    assert deletion_matrix.iloc[3, 8] == 1
    assert deletion_matrix.iloc[4, 9] == 1
    
    # empty value, villed with 0
    #column 0
    assert deletion_matrix.iloc[0, 0] == 0
    assert deletion_matrix.iloc[1, 0] == 0
    assert deletion_matrix.iloc[2, 0] == 0
    assert deletion_matrix.iloc[3, 0] == 0
    
    #column 1
    assert deletion_matrix.iloc[0, 1] == 0
    assert deletion_matrix.iloc[1, 1] == 0
    assert deletion_matrix.iloc[2, 1] == 0
    assert deletion_matrix.iloc[4, 1] == 0
    
    #columns 2
    assert deletion_matrix.iloc[0, 2] == 0
    assert deletion_matrix.iloc[1, 2] == 0
    assert deletion_matrix.iloc[3, 2] == 0
    assert deletion_matrix.iloc[4, 2] == 0
    
    #columns 3
    assert deletion_matrix.iloc[0, 3] == 0
    assert deletion_matrix.iloc[2, 3] == 0
    assert deletion_matrix.iloc[3, 3] == 0
    assert deletion_matrix.iloc[4, 3] == 0
    
    #columns 4
    assert deletion_matrix.iloc[1, 4] == 0
    assert deletion_matrix.iloc[2, 4] == 0
    assert deletion_matrix.iloc[3, 4] == 0
    assert deletion_matrix.iloc[4, 4] == 0
    
    #columns 5
    assert deletion_matrix.iloc[1, 5] == 0
    assert deletion_matrix.iloc[2, 5] == 0
    assert deletion_matrix.iloc[3, 5] == 0
    assert deletion_matrix.iloc[4, 5] == 0
    
    #columns 6
    assert deletion_matrix.iloc[0, 6] == 0
    assert deletion_matrix.iloc[2, 6] == 0
    assert deletion_matrix.iloc[3, 6] == 0
    assert deletion_matrix.iloc[4, 6] == 0
    
    #columns 7
    assert deletion_matrix.iloc[0, 7] == 0
    assert deletion_matrix.iloc[1, 7] == 0
    assert deletion_matrix.iloc[3, 7] == 0
    assert deletion_matrix.iloc[4, 7] == 0
    
    #columns 8
    assert deletion_matrix.iloc[0, 8] == 0
    assert deletion_matrix.iloc[1, 8] == 0
    assert deletion_matrix.iloc[2, 8] == 0
    assert deletion_matrix.iloc[4, 8] == 0
    
    #columns 9
    assert deletion_matrix.iloc[0, 9] == 0
    assert deletion_matrix.iloc[1, 9] == 0
    assert deletion_matrix.iloc[2, 9] == 0
    assert deletion_matrix.iloc[3, 9] == 0
    
    return deletion_matrix
    




def test_create_plots_create_matrix_deletion_percent():
    #ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("15_dels_ladder_for_matrix.fasta")
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("15_dels_ladder_for_matrix.fasta")
    cov = count_indels._get_coverage("15_dels_ladder_for_matrix.fasta")
    total_dels = count_indels._count_dels("15_dels_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("15_dels_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    

    df_dels = count_indels._create_df(ref_seq, total_dels)
    df_ins = count_indels._create_df(ref_seq, total_ins)
    df_cov = count_indels._create_df_cov(cov)
    
    count_indels._save_df(df_dels, df_ins, df_cov, "15_dels_ladder_for_matrix.fasta")
        
    deletion_matrix = create_matrices._create_matrix("15_dels_ladder_for_matrix.xlsx", "deletions")
    
    perc_deletion_matrix = create_matrices._create_matrix_percent("15_dels_ladder_for_matrix.xlsx",
                                                                      deletion_matrix,
                                                                      cov,
                                                                     "deletions")
    
    assert perc_deletion_matrix.iloc[4, 0] == 10
    assert perc_deletion_matrix.iloc[3, 1] == 10
    assert perc_deletion_matrix.iloc[2, 2] == 10
    assert perc_deletion_matrix.iloc[1, 3] == 10
    assert perc_deletion_matrix.iloc[0, 4] == 10
    assert perc_deletion_matrix.iloc[0, 5] == 10
    assert perc_deletion_matrix.iloc[1, 6] == 5
    assert perc_deletion_matrix.iloc[2, 7] == 5
    assert perc_deletion_matrix.iloc[3, 8] == 5
    assert perc_deletion_matrix.iloc[4, 9] == 5
    
    
    # empty value, villed with 0
    #column 0
    assert perc_deletion_matrix.iloc[0, 0] == 0
    assert perc_deletion_matrix.iloc[1, 0] == 0
    assert perc_deletion_matrix.iloc[2, 0] == 0
    assert perc_deletion_matrix.iloc[3, 0] == 0
    
    #column 1
    assert perc_deletion_matrix.iloc[0, 1] == 0
    assert perc_deletion_matrix.iloc[1, 1] == 0
    assert perc_deletion_matrix.iloc[2, 1] == 0
    assert perc_deletion_matrix.iloc[4, 1] == 0
    
    #columns 2
    assert perc_deletion_matrix.iloc[0, 2] == 0
    assert perc_deletion_matrix.iloc[1, 2] == 0
    assert perc_deletion_matrix.iloc[3, 2] == 0
    assert perc_deletion_matrix.iloc[4, 2] == 0
    
    #columns 3
    assert perc_deletion_matrix.iloc[0, 3] == 0
    assert perc_deletion_matrix.iloc[2, 3] == 0
    assert perc_deletion_matrix.iloc[3, 3] == 0
    assert perc_deletion_matrix.iloc[4, 3] == 0
    
    #columns 4
    assert perc_deletion_matrix.iloc[1, 4] == 0
    assert perc_deletion_matrix.iloc[2, 4] == 0
    assert perc_deletion_matrix.iloc[3, 4] == 0
    assert perc_deletion_matrix.iloc[4, 4] == 0
    
    #columns 5
    assert perc_deletion_matrix.iloc[1, 5] == 0
    assert perc_deletion_matrix.iloc[2, 5] == 0
    assert perc_deletion_matrix.iloc[3, 5] == 0
    assert perc_deletion_matrix.iloc[4, 5] == 0
    
    #columns 6
    assert perc_deletion_matrix.iloc[0, 6] == 0
    assert perc_deletion_matrix.iloc[2, 6] == 0
    assert perc_deletion_matrix.iloc[3, 6] == 0
    assert perc_deletion_matrix.iloc[4, 6] == 0
    
    #columns 7
    assert perc_deletion_matrix.iloc[0, 7] == 0
    assert perc_deletion_matrix.iloc[1, 7] == 0
    assert perc_deletion_matrix.iloc[3, 7] == 0
    assert perc_deletion_matrix.iloc[4, 7] == 0
    
    #columns 8
    assert perc_deletion_matrix.iloc[0, 8] == 0
    assert perc_deletion_matrix.iloc[1, 8] == 0
    assert perc_deletion_matrix.iloc[2, 8] == 0
    assert perc_deletion_matrix.iloc[4, 8] == 0
    
    #columns 9
    assert perc_deletion_matrix.iloc[0, 9] == 0
    assert perc_deletion_matrix.iloc[1, 9] == 0
    assert perc_deletion_matrix.iloc[2, 9] == 0
    assert perc_deletion_matrix.iloc[3, 9] == 0
    
    return perc_deletion_matrix
    

    

def test_create_plots_create_matrix_insertions_raw_count():
    """
    """
    #ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("16_ins_ladder_for_matrix.fasta")
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("16_ins_ladder_for_matrix.fasta")
    cov = count_indels._get_coverage("16_ins_ladder_for_matrix.fasta")
    total_dels = count_indels._count_dels("16_ins_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("16_ins_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    
    
    df_dels = count_indels._create_df(ref_seq, total_dels)
    df_ins = count_indels._create_df(ref_seq, total_ins)
    df_cov = count_indels._create_df_cov(cov)
    
    count_indels._save_df(df_dels, df_ins, df_cov, "16_ins_ladder_for_matrix.fasta")
    
    insertion_matrix = create_matrices._create_matrix("16_ins_ladder_for_matrix.xlsx", "insertions")
    try:
        deletions_matrix = create_matrices._create_matrix("16_ins_ladder_for_matrix.xlsx", "insertions")
    except ValueError:
        print("no insertion matrix")
        
    assert insertion_matrix.iloc[3, 1] == 3
    assert insertion_matrix.iloc[2, 2] == 2
    assert insertion_matrix.iloc[1, 3] == 2
    assert insertion_matrix.iloc[0, 4] == 2
    assert insertion_matrix.iloc[1, 5] == 1
    
    
    # column 0
    assert insertion_matrix.iloc[0, 0] == 0
    assert insertion_matrix.iloc[1, 0] == 0
    assert insertion_matrix.iloc[2, 0] == 0
    assert insertion_matrix.iloc[3, 0] == 0
    
    # column 1
    assert insertion_matrix.iloc[0, 1] == 0
    assert insertion_matrix.iloc[1, 1] == 0
    assert insertion_matrix.iloc[2, 1] == 0
    
    # column 2
    assert insertion_matrix.iloc[0, 2] == 0
    assert insertion_matrix.iloc[1, 2] == 0
    assert insertion_matrix.iloc[3, 2] == 0
   
    # column 3
    assert insertion_matrix.iloc[0, 3] == 0
    assert insertion_matrix.iloc[2, 3] == 0
    assert insertion_matrix.iloc[3, 3] == 0
    
    # column 4
    assert insertion_matrix.iloc[1, 4] == 0
    assert insertion_matrix.iloc[2, 4] == 0
    assert insertion_matrix.iloc[3, 4] == 0
    
    
    return insertion_matrix
    




def test_create_plots_create_matrix_insertions_percent():
    #ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("16_ins_ladder_for_matrix.fasta")
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("16_ins_ladder_for_matrix.fasta")
    cov = count_indels._get_coverage("16_ins_ladder_for_matrix.fasta")
    total_dels = count_indels._count_dels("16_ins_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("16_ins_ladder_for_matrix.fasta", ref_seq, ref_seq_id, cov)
    
    df_dels = count_indels._create_df(ref_seq, total_dels)
    df_ins = count_indels._create_df(ref_seq, total_ins)
    df_cov = count_indels._create_df_cov(cov)
    
    count_indels._save_df(df_dels, df_ins, df_cov, "16_ins_ladder_for_matrix.fasta")
    
    insertion_matrix = create_matrices._create_matrix("16_ins_ladder_for_matrix.xlsx", "insertions")
    
    insertion_matrix_percent = create_matrices._create_matrix_percent("16_ins_ladder_for_matrix.xlsx",
                                                                      insertion_matrix,
                                                                      cov,
                                                                     "insertions")
   
    assert insertion_matrix_percent.iloc[3, 1] == 30
    assert insertion_matrix_percent.iloc[2, 2] == 20
    assert insertion_matrix_percent.iloc[1, 3] == 20
    assert insertion_matrix_percent.iloc[0, 4] == 20
    assert insertion_matrix_percent.iloc[1, 5] == 10
    
    
    
    # column 0
    assert insertion_matrix_percent.iloc[0, 0] == 0
    assert insertion_matrix_percent.iloc[1, 0] == 0
    assert insertion_matrix_percent.iloc[2, 0] == 0
    assert insertion_matrix_percent.iloc[3, 0] == 0
    
    # column 1
    assert insertion_matrix_percent.iloc[0, 1] == 0
    assert insertion_matrix_percent.iloc[1, 1] == 0
    assert insertion_matrix_percent.iloc[2, 1] == 0
    
    # column 2
    assert insertion_matrix_percent.iloc[0, 2] == 0
    assert insertion_matrix_percent.iloc[1, 2] == 0
    assert insertion_matrix_percent.iloc[3, 2] == 0
   
    # column 3
    assert insertion_matrix_percent.iloc[0, 3] == 0
    assert insertion_matrix_percent.iloc[2, 3] == 0
    assert insertion_matrix_percent.iloc[3, 3] == 0
    
    # column 4
    assert insertion_matrix_percent.iloc[1, 4] == 0
    assert insertion_matrix_percent.iloc[2, 4] == 0
    assert insertion_matrix_percent.iloc[3, 4] == 0
    
    return insertion_matrix_percent
    


# calls to have a look at the data
#insertion_matrix = test_create_plots_create_matrix_insertions_raw_count()
#print("insertion matrix: ", insertion_matrix, sep="\n")


#insertion_matrix_percent = test_create_plots_create_matrix_insertions_percent()
#print("insertion matrix percent: ", insertion_matrix_percent, sep="\n")


#deletion_matrix = test_create_plots_create_matrix_deletion_raw_count()
#print("deletion matrix: ", deletion_matrix)


#deletions_matrix_percent = test_create_plots_create_matrix_deletion_percent()
#print("deletion matrix percent: ", deletions_matrix_percent)









