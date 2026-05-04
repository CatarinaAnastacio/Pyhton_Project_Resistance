import pytest

from resistance import reverse_com_strand


@pytest.mark.parametrize("seq,expected", [("ATCG", "CGAT"),
    ("AAAA", "TTTT"),("GCGC", "GCGC"),("A", "T"),("", "")])

def test_reverse_com_strand_normal(seq, expected):
    assert reverse_com_strand(seq) == expected


# Test to confirm that N continues present in the reverse complement
def test_reverse_com_strand_with_N():
    assert reverse_com_strand("ATCGN") == "NCGAT"

# Test to confirm that characters that are not ACGT become N
def test_reverse_com_strand_non_base():
    assert reverse_com_strand("ATCGX") == "NCGAT"

# Test to check if Reverse complement appplied twice returns the same sequence
def test_double_reverse_com_strand():
    seq = "ATCGATCGATCGATCG"
    assert reverse_com_strand(reverse_com_strand(seq)) == seq



###### STILL BASIC - DO WE NEED MORE????