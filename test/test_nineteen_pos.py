
import pytest

from resistance import nineteen_pos

# Test to confirm that seqeunces shorter than 19 should yield nothing
def test_nineteen_pos_short_seq():
    result = list(nineteen_pos("ATCG"))
    assert result == []


# Test to confirm that a sequence of length 19 yields only one k-mer in position 0
def test_nineteen_pos_length_19():
    seq = "ATGCGTACGTTAGCCTGAA"
    result = list(nineteen_pos(seq))
    assert result == [(0, "ATGCGTACGTTAGCCTGAA")]


# Test to confirm that a sequence of length 20 yields two k-mer in position 0 and 1
def test_nineteen_pos_length_20():
    seq = "CAAAAAAAAATTTTTTTTTC"
    result = list(nineteen_pos(seq))
    assert len(result) == 2
    assert result[0] == (0, "CAAAAAAAAATTTTTTTTT")
    assert result[1] == (1, "AAAAAAAAATTTTTTTTTC")


# Test to confirm that k-mers with N are skipped
def test_nineteen_pos_with_N():
    seq = "AAAAAAAAAAAAAAAAAAA" + "N" + "AAAAAAAAAAAAAAAAAAA" 
    # Positions 0 and 20 are valid, positions 1 until 19 contain the N
    result = list(nineteen_pos(seq))
    positions = [pos for pos, _ in result]
    assert positions == [0, 20]
    # There should be no k-mer with N
    for _, kmer in result:
        assert "N" not in kmer


# Test to confirm the number of k-mers for a normal sequence without N
def test_nineteen_pos_count():
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # length 50
    result = list(nineteen_pos(seq))
    assert len(result) == 32