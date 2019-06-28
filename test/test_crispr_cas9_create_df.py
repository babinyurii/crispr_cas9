# -*- coding: utf-8 -*-
"""
Created on Tue May 21 18:03:12 2019

@author: babin
"""

import sys
sys.path.append("..")

from crispr_cas9 import count_indels


def test_count_indels_single_del_from_head():
    """single deletions over the whole ref_seq
    no insertions"""
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("13_df_dels.fasta")
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("13_df_dels.fasta")
    cov = count_indels._get_coverage("13_df_dels.fasta")
    total_dels = count_indels._count_dels("13_df_dels.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("13_df_dels.fasta", ref_seq, ref_seq_id, cov)
    
    df_dels = count_indels._create_df(ref_seq, total_dels)
    df_ins = count_indels._create_df(ref_seq, total_ins)
    df_cov = count_indels._create_df_cov(cov)
    
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
    
    assert df_cov["coverage"].item() == 1

    
    return df_dels, df_ins, df_cov



def test_count_indels_insertion_ladder_forward():
    """
    """
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("10_ins_ladder_forward.fasta")
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("10_ins_ladder_forward.fasta")
    cov = count_indels._get_coverage("10_ins_ladder_forward.fasta")
    total_dels = count_indels._count_dels("10_ins_ladder_forward.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("10_ins_ladder_forward.fasta", ref_seq, ref_seq_id, cov)
    
    df_dels = count_indels._create_df(ref_seq, total_dels)
    df_ins = count_indels._create_df(ref_seq, total_ins)
    df_cov = count_indels._create_df_cov(cov)
    
    
    assert df_ins[("A", 1)].isna().item() == True
    assert df_ins[("T", 2)].item() == 1
    assert df_ins[("G", 3)].item() == 2
    assert df_ins[("C", 4)].item() == 3
    assert df_ins[("A", 5)].item() == 4
    assert df_ins[("T", 6)].item() == 3
    
    assert df_dels.empty == True
    
    assert df_cov["coverage"].item() == 5
    
    
    return df_dels, df_ins, df_cov





#df_dels, df_ins, df_cov = test_count_indels_insertion_ladder_forward()

#print(df_ins)
#print(df_cov)


#df_dels, df_ins, df_cov  = test_count_indels_single_del_from_head()

#print(df_dels)
#print(df_ins)
#print(df_cov)