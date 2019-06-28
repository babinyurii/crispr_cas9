# -*- coding: utf-8 -*-
"""
Created on Mon May 20 10:34:12 2019

@author: babin
"""
import sys
sys.path.append("..")

from crispr_cas9 import count_indels
 
"""
nested arrays total_dels, total_ins have the following structure:
- each array entry corresponds to the position according to
the reference sequence
- integers in the entry are the single indels events, and the integer is the indel length
f.e., if we have total_dels like: [[1, 2, 3], [...], ...],  it means that at the first 
position there are there deletions of 1, 2 and 3 nucleotide lengths
"""    
def test_count_indels_single_del_from_head():
    """single deletions over the whole ref_seq
    no insertions"""
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("1_dels_single_del_from_head.fasta")
    cov = count_indels._get_coverage("1_dels_single_del_from_head.fasta")
    total_dels = count_indels._count_dels("1_dels_single_del_from_head.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("1_dels_single_del_from_head.fasta", ref_seq, ref_seq_id, cov)
    
    assert len(ref_seq) == 10
    assert cov == 10
    
    assert total_dels[0][0] == 1
    assert total_dels[1][0] == 1
    assert total_dels[2][0] == 1
    assert total_dels[3][0] == 1
    assert total_dels[4][0] == 1
    assert total_dels[5][0] == 1
    assert total_dels[6][0] == 1
    assert total_dels[7][0] == 1
    assert total_dels[8][0] == 1
    assert total_dels[9][0] == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    assert total_ins[3] == []
    assert total_ins[4] == []
    assert total_ins[5] == []
    assert total_ins[6] == []
    assert total_ins[7] == []
    assert total_ins[8] == []
    assert total_ins[9] == []
    
    
    
def test_count_indels_single_del_from_tail():

    ref_seq, ref_seq_id = count_indels._get_ref_seq("2_dels_single_del_from_tail.fasta")
    cov = count_indels._get_coverage("2_dels_single_del_from_tail.fasta")
    total_dels = count_indels._count_dels("2_dels_single_del_from_tail.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("2_dels_single_del_from_tail.fasta", ref_seq, ref_seq_id, cov)


    assert total_dels[0][0] == 1
    assert total_dels[1][0] == 1
    assert total_dels[2][0] == 1
    assert total_dels[3][0] == 1
    assert total_dels[4][0] == 1
    assert total_dels[5][0] == 1
    assert total_dels[6][0] == 1
    assert total_dels[7][0] == 1
    assert total_dels[8][0] == 1
    assert total_dels[9][0] == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    assert total_ins[3] == []
    assert total_ins[4] == []
    assert total_ins[5] == []
    assert total_ins[6] == []
    assert total_ins[7] == []
    assert total_ins[8] == []
    assert total_ins[9] == []
    
 
    
def test_count_indels_single_del_hang_at_head():

    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("3_dels_single_del_hang_at_head.fasta")
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("3_dels_single_del_hang_at_head.fasta")
    cov = count_indels._get_coverage("3_dels_single_del_hang_at_head.fasta")
    total_dels = count_indels._count_dels("3_dels_single_del_hang_at_head.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("3_dels_single_del_hang_at_head.fasta", ref_seq, ref_seq_id, cov)
    
    assert total_dels[0][0] == 1
    assert total_dels[0][1] == 1
    assert total_dels[0][2] == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    


def test_count_indels_single_del_hang_at_tail():

    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("4_dels_single_del_hang_at_tail.fasta")
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("4_dels_single_del_hang_at_tail.fasta")
    cov = count_indels._get_coverage("4_dels_single_del_hang_at_tail.fasta")
    total_dels = count_indels._count_dels("4_dels_single_del_hang_at_tail.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("4_dels_single_del_hang_at_tail.fasta", ref_seq, ref_seq_id, cov)
    
    
    
    assert total_dels[9][0] == 1
    assert total_dels[9][1] == 1
    assert total_dels[9][2] == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    



def test_count_indels_hang_from_head_various_len():
 
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("5_dels_hanging_from_head_var_len.fasta")
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("5_dels_hanging_from_head_var_len.fasta")
    cov = count_indels._get_coverage("5_dels_hanging_from_head_var_len.fasta")
    total_dels = count_indels._count_dels("5_dels_hanging_from_head_var_len.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("5_dels_hanging_from_head_var_len.fasta", ref_seq, ref_seq_id, cov)
    
    
    
    assert total_dels[0][0] == 1
    assert total_dels[0][1] == 2
    assert total_dels[0][2] == 3
    assert total_dels[0][3] == 4
    assert total_dels[0][4] == 5
    assert total_dels[0][5] == 6
    assert total_dels[0][6] == 7
    assert total_dels[0][7] == 8
    assert total_dels[0][8] == 9
    
    # all the other entries are empty
    assert total_dels[1] == []
    assert total_dels[2] == []
    assert total_dels[3] == []
    assert total_dels[4] == []
    assert total_dels[5] == []
    assert total_dels[6] == []
    assert total_dels[7] == []
    assert total_dels[8] == []
    assert total_dels[9] == []
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    assert total_ins[3] == []
    assert total_ins[4] == []
    assert total_ins[5] == []
    assert total_ins[6] == []
    assert total_ins[7] == []
    assert total_ins[8] == []
    assert total_ins[9] == []
    

def test_count_indels_hang_from_tail_various_len():
  
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("6_dels_hanging_from_tail_var_len.fasta")
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("6_dels_hanging_from_tail_var_len.fasta")
    cov = count_indels._get_coverage("6_dels_hanging_from_tail_var_len.fasta")
    total_dels = count_indels._count_dels("6_dels_hanging_from_tail_var_len.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("6_dels_hanging_from_tail_var_len.fasta", ref_seq, ref_seq_id, cov)
    
    
    #assert total_dels[0][0] == 10
    assert total_dels[1][0] == 9
    assert total_dels[2][0] == 8
    assert total_dels[3][0] == 7
    assert total_dels[4][0] == 6
    assert total_dels[5][0] == 5
    assert total_dels[6][0] == 4
    assert total_dels[7][0] == 3
    assert total_dels[8][0] == 2
    assert total_dels[9][0] == 1
    
    # only one del at each entry
    assert len(total_dels[0]) == 0
    assert len(total_dels[1]) == 1
    assert len(total_dels[2]) == 1
    assert len(total_dels[3]) == 1
    assert len(total_dels[4]) == 1
    assert len(total_dels[5]) == 1
    assert len(total_dels[6]) == 1
    assert len(total_dels[7]) == 1
    assert len(total_dels[8]) == 1
    assert len(total_dels[9]) == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    assert total_ins[3] == []
    assert total_ins[4] == []
    assert total_ins[5] == []
    assert total_ins[6] == []
    assert total_ins[7] == []
    assert total_ins[8] == []
    assert total_ins[9] == []
    



def test_count_indels_del_ladder():
  
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("7_dels_ladder.fasta")
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("7_dels_ladder.fasta")
    cov = count_indels._get_coverage("7_dels_ladder.fasta")
    total_dels = count_indels._count_dels("7_dels_ladder.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("7_dels_ladder.fasta", ref_seq, ref_seq_id, cov)
    
    
    assert total_dels[0][0] == 1
    assert total_dels[1][0] == 2
    assert total_dels[2][0] == 3
    assert total_dels[3][0] == 4
    assert total_dels[4][0] == 5
    assert total_dels[5][0] == 5
    assert total_dels[6][0] == 4
    assert total_dels[7][0] == 3
    assert total_dels[8][0] == 2
    assert total_dels[9][0] == 1
    
    # only one del at each entry
    assert len(total_dels[0]) == 1
    assert len(total_dels[1]) == 1
    assert len(total_dels[2]) == 1
    assert len(total_dels[3]) == 1
    assert len(total_dels[4]) == 1
    assert len(total_dels[5]) == 1
    assert len(total_dels[6]) == 1
    assert len(total_dels[7]) == 1
    assert len(total_dels[8]) == 1
    assert len(total_dels[9]) == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    assert total_ins[3] == []
    assert total_ins[4] == []
    assert total_ins[5] == []
    assert total_ins[6] == []
    assert total_ins[7] == []
    assert total_ins[8] == []
    assert total_ins[9] == []
    
   

    
"""
insertions testing
"""
  
def test_count_indels_single_insertions_from_head():
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("8_ins_single_ins_from_head.fasta")
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("8_ins_single_ins_from_head.fasta")
    cov = count_indels._get_coverage("8_ins_single_ins_from_head.fasta")
    total_dels = count_indels._count_dels("8_ins_single_ins_from_head.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("8_ins_single_ins_from_head.fasta", ref_seq, ref_seq_id, cov)
    
    
    
    assert total_ins[1][0] == 1
    assert total_ins[2][0] == 1
    assert total_ins[3][0] == 1
    assert total_ins[4][0] == 1
    assert total_ins[5][0] == 1
    assert total_ins[6][0] == 1
    assert total_ins[7][0] == 1
    assert total_ins[8][0] == 1
    assert total_ins[9][0] == 1
    
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 1
    assert len(total_ins[2]) == 1
    assert len(total_ins[3]) == 1
    assert len(total_ins[4]) == 1
    assert len(total_ins[5]) == 1
    assert len(total_ins[6]) == 1
    assert len(total_ins[7]) == 1
    assert len(total_ins[8]) == 1
    assert len(total_ins[9]) == 1
   
    
    
    assert total_dels[0] == []
    assert total_dels[1] == []
    assert total_dels[2] == []
    assert total_dels[3] == []
    assert total_dels[4] == []
    assert total_dels[5] == []
    assert total_dels[6] == []
    assert total_dels[7] == []
    assert total_dels[8] == []
    assert total_dels[9] == []
    
    
def test_count_indels_single_insertions_from_tail():
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("9_ins_single_ins_from_tail.fasta")
    
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("9_ins_single_ins_from_tail.fasta")
    cov = count_indels._get_coverage("9_ins_single_ins_from_tail.fasta")
    total_dels = count_indels._count_dels("9_ins_single_ins_from_tail.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("9_ins_single_ins_from_tail.fasta", ref_seq, ref_seq_id, cov)
    
    assert total_ins[1][0] == 1
    assert total_ins[2][0] == 1
    assert total_ins[3][0] == 1
    assert total_ins[4][0] == 1
    assert total_ins[5][0] == 1
    assert total_ins[6][0] == 1
    assert total_ins[7][0] == 1
    assert total_ins[8][0] == 1
    assert total_ins[9][0] == 1
    
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 1
    assert len(total_ins[2]) == 1
    assert len(total_ins[3]) == 1
    assert len(total_ins[4]) == 1
    assert len(total_ins[5]) == 1
    assert len(total_ins[6]) == 1
    assert len(total_ins[7]) == 1
    assert len(total_ins[8]) == 1
    assert len(total_ins[9]) == 1
   
    assert total_dels[0] == []
    assert total_dels[1] == []
    assert total_dels[2] == []
    assert total_dels[3] == []
    assert total_dels[4] == []
    assert total_dels[5] == []
    assert total_dels[6] == []
    assert total_dels[7] == []
    assert total_dels[8] == []
    assert total_dels[9] == []    
    
    

def test_count_indels_insertions_var_lens_ladder_forward():
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("10_ins_ladder_forward.fasta")
    
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("10_ins_ladder_forward.fasta")
    cov = count_indels._get_coverage("10_ins_ladder_forward.fasta")
    total_dels = count_indels._count_dels("10_ins_ladder_forward.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("10_ins_ladder_forward.fasta", ref_seq, ref_seq_id, cov)
    
    #assert total_ins[0][0] == 0
    assert total_ins[1][0] == 1
    assert total_ins[2][0] == 2
    assert total_ins[3][0] == 3
    assert total_ins[4][0] == 4
    assert total_ins[5][0] == 3
    
    # check if there's one insertion under each position
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 1
    assert len(total_ins[2]) == 1
    assert len(total_ins[3]) == 1
    assert len(total_ins[4]) == 1
    assert len(total_ins[5]) == 1
    

    assert total_dels[0] == []
    assert total_dels[1] == []
    assert total_dels[2] == []
    assert total_dels[3] == []
    assert total_dels[4] == []
    assert total_dels[5] == []
      
    
    
    


def test_count_indels_insertions_var_lens_ladder():
    #ref_seq, total_dels, total_ins, cov = count_indels._count_indels("11_ins_ladder.fasta")
    
    ref_seq, ref_seq_id = count_indels._get_ref_seq("11_ins_ladder.fasta")
    cov = count_indels._get_coverage("11_ins_ladder.fasta")
    total_dels = count_indels._count_dels("11_ins_ladder.fasta", ref_seq, ref_seq_id, cov)
    total_ins = count_indels._count_ins("11_ins_ladder.fasta", ref_seq, ref_seq_id, cov)
    
    assert total_ins[1][0] == 1
    assert total_ins[1][1] == 1
    assert total_ins[2][0] == 2
    assert total_ins[2][1] == 2
    assert total_ins[3][0] == 3
    assert total_ins[3][1] == 3
    assert total_ins[4][0] == 4
    assert total_ins[4][1] == 4
    assert total_ins[5][0] == 3
    
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 2
    assert len(total_ins[2]) == 2
    assert len(total_ins[3]) == 2
    assert len(total_ins[4]) == 2
    assert len(total_ins[5]) == 1
    
    
    assert total_dels[0] == []
    assert total_dels[1] == []
    assert total_dels[2] == []
    assert total_dels[3] == []
    assert total_dels[4] == []
    assert total_dels[5] == []



