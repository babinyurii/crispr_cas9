# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:03:12 2019

@author: babin
"""
import numpy as np
import sys
sys.path.append("..")

from crispr_cas9 import crispr_count_indels


def test_count_indels_single_del_from_head():
    """single deletions over the whole ref_seq
    no insertions"""
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("13_df_dels.fasta")
    
    df_dels, df_ins, df_cov = crispr_count_indels._create_df(ref_seq,
                                                             total_dels,
                                                             total_ins, 
                                                             "13_df_dels.fasta",
                                                             cov)
    # we should use item() method to return Python scalar value
    # otherwise we'll get ValueError: The truth value of a Series is ambiguous.
    # to check up compare: 
    # type(df_dels[("A", 1)].isna()) and type(df_dels[("A", 1)].isna().item())
    # the first we'll be pandas.core.series.Series, the second - bool
    
    assert df_dels[("T", 2)].item() == 1
    assert df_dels[("A", 5)].item() == 2
 
    assert df_dels[("A", 1)].isna().item() == True
    assert df_dels[("G", 3)].isna().item() == True
    assert df_dels[("C", 4)].isna().item() == True
    assert df_dels[("T", 6)].isna().item() == True
    assert df_dels[("G", 7)].isna().item() == True

    assert df_ins.empty == True
    
    
    return df_dels, df_ins, df_cov


def test_count_indels_insertion_ladder_forward():
    """
    """
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("10_ins_ladder_forward.fasta")
    
    df_dels, df_ins, df_cov = crispr_count_indels._create_df(ref_seq,
                                                             total_dels,
                                                             total_ins, 
                                                             "10_ins_ladder_forward.fasta",
                                                             cov)
    assert df_ins[("A", 1)].isna().item() == True
    assert df_ins[("T", 2)].item() == 1
    assert df_ins[("G", 3)].item() == 2
    assert df_ins[("C", 4)].item() == 3
    assert df_ins[("A", 5)].item() == 4
    assert df_ins[("T", 6)].item() == 3
    
    assert df_dels.empty == True
    
    return df_ins


df_ins = test_count_indels_insertion_ladder_forward()

print(df_ins)


#df_dels, df_ins, df_cov  = test_count_indels_single_del_from_head()

#print(df_dels)
#print(df_ins)
#print(df_cov)