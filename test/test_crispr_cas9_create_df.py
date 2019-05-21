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
    
    assert int(df_dels[("T", 2)].item()) == 1
    assert int(df_dels[("A", 5)].item()) == 2
    
    #assert df_dels[("A", 1)].item().isna() == True
    
    
    
    return df_dels, df_ins, df_cov


df_dels, df_ins, df_cov  = test_count_indels_single_del_from_head()

print(df_dels)
print(df_ins)
print(df_cov)