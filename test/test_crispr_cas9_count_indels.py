# -*- coding: utf-8 -*-
"""
Created on Mon May 20 10:34:12 2019

@author: babin
"""
import sys
sys.path.append("..")

from crispr_cas9 import crispr_count_indels
 
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
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("1_dels_single_del_from_head.fasta")
    
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
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("2_dels_single_del_from_tail.fasta")

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

    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("3_dels_single_del_hang_at_head.fasta")
    
    assert total_dels[0][0] == 1
    assert total_dels[0][1] == 1
    assert total_dels[0][2] == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    


def test_count_indels_single_del_hang_at_tail():

    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("4_dels_single_del_hang_at_tail.fasta")
    
    assert total_dels[9][0] == 1
    assert total_dels[9][1] == 1
    assert total_dels[9][2] == 1
    
    assert total_ins[0] == []
    assert total_ins[1] == []
    assert total_ins[2] == []
    



def test_count_indels_hang_from_head_various_len():
 
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("5_dels_hanging_from_head_var_len.fasta")
    
    assert total_dels[0][0] == 1
    assert total_dels[0][1] == 2
    assert total_dels[0][2] == 3
    assert total_dels[0][3] == 4
    assert total_dels[0][4] == 5
    assert total_dels[0][5] == 6
    assert total_dels[0][6] == 7
    assert total_dels[0][7] == 8
    assert total_dels[0][8] == 9
    assert total_dels[0][9] == 10
    
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
  
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("6_dels_hanging_from_tail_var_len.fasta")
    
    assert total_dels[0][0] == 10
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
    



def test_count_indels_del_ladder():
  
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("7_dels_ladder.fasta")
    
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
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("8_ins_single_ins_from_head.fasta")
    
    assert total_ins[1][0] == 1
    assert total_ins[3][0] == 1
    assert total_ins[5][0] == 1
    assert total_ins[7][0] == 1
    assert total_ins[9][0] == 1
    assert total_ins[11][0] == 1
    assert total_ins[13][0] == 1
    assert total_ins[15][0] == 1
    assert total_ins[17][0] == 1
    
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 1
    assert len(total_ins[2]) == 0
    assert len(total_ins[3]) == 1
    assert len(total_ins[4]) == 0
    assert len(total_ins[5]) == 1
    assert len(total_ins[6]) == 0
    assert len(total_ins[7]) == 1
    assert len(total_ins[8]) == 0
    assert len(total_ins[9]) == 1
    assert len(total_ins[10]) == 0
    assert len(total_ins[11]) == 1
    assert len(total_ins[12]) == 0
    assert len(total_ins[13]) == 1
    assert len(total_ins[14]) == 0
    assert len(total_ins[15]) == 1
    assert len(total_ins[16]) == 0
    assert len(total_ins[17]) == 1
    
    
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
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("9_ins_single_ins_from_tail.fasta")
    
    assert total_ins[1][0] == 1
    assert total_ins[3][0] == 1
    assert total_ins[5][0] == 1
    assert total_ins[7][0] == 1
    assert total_ins[9][0] == 1
    assert total_ins[11][0] == 1
    assert total_ins[13][0] == 1
    assert total_ins[15][0] == 1
    assert total_ins[17][0] == 1
    
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 1
    assert len(total_ins[2]) == 0
    assert len(total_ins[3]) == 1
    assert len(total_ins[4]) == 0
    assert len(total_ins[5]) == 1
    assert len(total_ins[6]) == 0
    assert len(total_ins[7]) == 1
    assert len(total_ins[8]) == 0
    assert len(total_ins[9]) == 1
    assert len(total_ins[10]) == 0
    assert len(total_ins[11]) == 1
    assert len(total_ins[12]) == 0
    assert len(total_ins[13]) == 1
    assert len(total_ins[14]) == 0
    assert len(total_ins[15]) == 1
    assert len(total_ins[16]) == 0
    assert len(total_ins[17]) == 1
    
    
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
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("10_ins_ladder_forward.fasta")
    
    assert total_ins[1][0] == 1
    assert total_ins[3][0] == 2
    assert total_ins[6][0] == 3
    assert total_ins[10][0] == 4
    assert total_ins[15][0] == 3
    
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 1
    assert len(total_ins[2]) == 0
    assert len(total_ins[3]) == 1
    assert len(total_ins[4]) == 0
    assert len(total_ins[5]) == 0
    assert len(total_ins[6]) == 1
    assert len(total_ins[7]) == 0
    assert len(total_ins[8]) == 0
    assert len(total_ins[9]) == 0
    assert len(total_ins[10]) == 1
    assert len(total_ins[11]) == 0
    assert len(total_ins[12]) == 0
    assert len(total_ins[13]) == 0
    assert len(total_ins[14]) == 0
    assert len(total_ins[15]) == 1
    assert len(total_ins[16]) == 0
    assert len(total_ins[17]) == 0
    assert len(total_ins[18]) == 0

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
    
    
    


def test_count_indels_insertions_var_lens_ladder():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("11_ins_ladder.fasta")
    
    assert total_ins[1][0] == 1
    assert total_ins[1][1] == 1
    assert total_ins[3][0] == 2
    assert total_ins[3][1] == 2
    assert total_ins[6][0] == 3
    assert total_ins[6][1] == 3
    assert total_ins[10][0] == 4
    assert total_ins[10][1] == 4
    assert total_ins[15][0] == 3
    
    assert len(total_ins[0]) == 0
    assert len(total_ins[1]) == 2
    assert len(total_ins[2]) == 0
    assert len(total_ins[3]) == 2
    assert len(total_ins[4]) == 0
    assert len(total_ins[5]) == 0
    assert len(total_ins[6]) == 2
    assert len(total_ins[7]) == 0
    assert len(total_ins[8]) == 0
    assert len(total_ins[9]) == 0
    assert len(total_ins[10]) == 2
    assert len(total_ins[11]) == 0
    assert len(total_ins[12]) == 0
    assert len(total_ins[13]) == 0
    assert len(total_ins[14]) == 0
    assert len(total_ins[15]) == 1
    assert len(total_ins[16]) == 0
    assert len(total_ins[17]) == 0
    assert len(total_ins[18]) == 0
    
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
    

def test_count_indels_dels_and_ins():
    ref_seq, total_dels, total_ins, cov  = crispr_count_indels._count_indels("12_dels_and_ins.fasta")
    
    assert total_dels[0][0] == 1
    assert total_dels[1] == []
    assert total_dels[2][0] == 1
    assert total_dels[3] == []
    assert total_dels[4] == []
    assert total_dels[5][0] == 1
    assert total_dels[6] == []
    assert total_dels[7] == []
    assert total_dels[8] == []
    assert total_dels[9][0] == 1
    
    assert total_ins[0] == []
    assert total_ins[1][0] == 1
    assert total_ins[2] == []
    assert total_ins[3][0] == 2
    assert total_ins[4] == []
    assert total_ins[5] == []
    assert total_ins[6][0] == 3
    assert total_ins[7] == []
    assert total_ins[8] == []
    assert total_ins[9] == []
    
    assert len(total_dels[0]) == 1
    assert len(total_dels[2]) == 1
    assert len(total_dels[5]) == 1
    assert len(total_dels[9]) == 1
    assert len(total_ins[1]) == 1
    assert len(total_ins[3]) == 1
    assert len(total_ins[6]) == 1








"""
just to have visual representation of the test data
"""

print("=========== single from head : dels_single_del_from_head.fasta =============")
ref_seq, total_dels, total_ins, cov  = crispr_count_indels._count_indels("1_dels_single_del_from_head.fasta")
print(total_dels)
print(total_ins)


print("=========== 12_dels_and_ins.fasta =============")
ref_seq, total_dels, total_ins, cov  = crispr_count_indels._count_indels("12_dels_and_ins.fasta")
print(total_dels)
print(total_ins)
















