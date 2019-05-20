# -*- coding: utf-8 -*-
"""
Created on Mon May 20 10:34:12 2019

@author: babin
"""
import sys
sys.path.append("..")

from crispr_cas9 import crispr_count_indels


def test_count_indels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("test.fasta")
    assert len(ref_seq) == 14
    
    assert total_ins[1][0] == 1, "position 1 has 1 ins of len 1"
    assert total_ins[3][0] == 1, "position 3 has 1 ins of len 1"
    assert total_dels[5][0] == 1, "position 5 has 1 del of len 1"
    
print("=========== single from head =============")
def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("dels_single_del_from_head.fasta")
    
    return total_dels, total_ins
    
    
t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)


print("============= single del from tail ============")
def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("dels_single_del_from_tail.fasta")
    
    return total_dels, total_ins

t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)



print("======= from head ==========")
def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("dels_hanging_from_head.fasta")
    
    return total_dels, total_ins
     
t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)

print("======= from tail ==========")

def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("dels_hanging_from_tail.fasta")
    
    return total_dels, total_ins
    
    
t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)




print("======= ladder ==========")

def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("dels_ladder.fasta")
    
    return total_dels, total_ins
    
    
t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)


print("======= insertions single ins from head ==========")

def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("ins_single_ins_from_head.fasta")
    
    return total_dels, total_ins
    
    
t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)

print("======= insertions single ins from tail ==========")

def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("ins_single_ins_from_tail.fasta")
    
    return total_dels, total_ins
    
    
t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)


print("======= insertions ladder ==========")

def test_count_indels_dels():
    ref_seq, total_dels, total_ins, cov = crispr_count_indels._count_indels("ins_ladder.fasta")
    
    return total_dels, total_ins
    
    
t_del, t_ins = test_count_indels_dels()
print(t_del)
print(t_ins)














