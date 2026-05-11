"""
Tests for the Fastq class.

The test data files contain the following:
    - kmer_match_testtdata.fastq.gz -> 5 reads, including a short read and a read with N's
    - empty_reads.fastq.gz -> empty file
    - no_match_reads.fastq.gz -> 3 reads that do not match any gene
Important: this tests use the gz folders instead of the simple .fastq files that the 'match' testes use
"""

import pytest

from src.fastq import Fastq


# Test to confirm the correct number of reads
def test_fastq_number_of_reads():

    fastq = Fastq("testdata/kmer_match_testtdata.fastq.gz")

    reads = list(fastq.fastqread())

    assert len(reads) == 5


# Test to confirm only DNA sequences are returned
def test_fastq_sequences_only():

    fastq = Fastq("testdata/kmer_match_testtdata.fastq.gz")

    reads = list(fastq.fastqread())

    for seq in reads:

        assert "@" not in seq
        assert "+" not in seq

        for base in seq:

            assert base in "ACGTN"


# Test to confirm the first read is correct
def test_fastq_first_read():

    fastq = Fastq("testdata/kmer_match_testtdata.fastq.gz")

    reads = list(fastq.fastqread())

    expected = "ACGTACGTACGTACGTACGTACGTACGTAC"

    assert reads[0] == expected


# Test to confirm empty FASTQ files return nothing
def test_fastq_empty_file():

    fastq = Fastq("testdata/empty_reads.fastq.gz")

    reads = list(fastq.fastqread())

    assert reads == []


# Test to confirm junk FASTQ reads are loaded
def test_fastq_junk_reads():

    fastq = Fastq("testdata/no_match_reads.fastq.gz")

    reads = list(fastq.fastqread())

    assert len(reads) == 3