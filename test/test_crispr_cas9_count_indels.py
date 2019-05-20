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
    