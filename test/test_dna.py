"""
Tests for the DNA class.
"""

import pytest

from src.dna import DNA

##### reverse_com_strand #####

# Test to confirm reverse complement of ACGT
def test_reverse_complement_basic():

    dna = DNA() 
    result = dna.reverse_com_strand("ACGT")
    
    assert result == "ACGT"


# Test to confirm reverse complement of AAAA
def test_reverse_complement_simple():
    
    dna = DNA()
    result = dna.reverse_com_strand("AAAA")
    
    assert result == "TTTT"


# Test to confirm reverse complement of ATCG
def test_reverse_complement_mixed():
    
    dna = DNA()
    result = dna.reverse_com_strand("ATCG")
    
    assert result == "CGAT"


# Test to confirm N bases become N
def test_reverse_complement_with_N():
    
    dna = DNA()
    result = dna.reverse_com_strand("ANCT")
    
    assert result == "AGNT"


# Test to confirm empty sequences return empty output
def test_reverse_complement_empty():
    
    dna = DNA()
    result = dna.reverse_com_strand("")
    
    assert result == ""


##### nineteen_pos #####

# Test to confirm sequences shorter than 19 return nothing
def test_nineteen_pos_short_sequence():
    
    dna = DNA()
    result = list(dna.nineteen_pos("ACGT"))
    
    assert result == []


# Test to confirm a sequence of length 19 returns one kmer
def test_nineteen_pos_exact_19():
    
    dna = DNA()
    seq = "AAAAAAAAAAAAAAAAAAA"
    
    result = list(dna.nineteen_pos(seq))

    assert len(result) == 1
    assert result[0] == (0, "AAAAAAAAAAAAAAAAAAA")


# Test to confirm correct number of kmers
def test_nineteen_pos_count():
    
    dna = DNA()
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAA"
    
    result = list(dna.nineteen_pos(seq))

    assert len(result) == 7


# Test to confirm kmers with N are skipped
def test_nineteen_pos_with_N():
    dna = DNA()
    seq = "AAAAAAAAAAAAAAAAAAA" + "N" + "AAAAAAAAAA"
    
    result = list(dna.nineteen_pos(seq))

    for _, kmer in result:
        assert "N" not in kmer


# Test to confirm positions increase correctly
def test_nineteen_pos_position_order():
    
    dna = DNA()
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    
    result = list(dna.nineteen_pos(seq))
    positions = [pos for pos, _ in result]

    assert positions == list(range(len(seq) - 18))