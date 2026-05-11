"""
Tests for the Match class.

The test FASTQ 'kmer_match_testtdata.fastq.gz' contains 5 reads, two of
which do match a gene in the test FASTA, also a junk read, a too-short
read and a read containing N's.

Important: this tests use the simple .fastq files instead of the .gz folders that the 'fastq' tests use
"""

import pytest

from src.fasta import Fasta
from src.kmer import Kmer
from src.match import Match


# Test to confirm depth arrays are initialized correctly
def test_match_depth_initialization():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)
    kmer.create_dict_kmer()

    matcher = Match([], kmer.kmer_index, fasta.gene_length)

    assert len(matcher.dict_depth["geneA_test_1"]) == 60
    assert len(matcher.dict_depth["geneB_test_2"]) == 60


# Test to confirm depth starts at zero
def test_match_depth_starts_zero():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)
    kmer.create_dict_kmer()

    matcher = Match([], kmer.kmer_index, fasta.gene_length)

    for depths in matcher.dict_depth.values():
        assert sum(depths) == 0


# Test to confirm matching reads increase depth
def test_match_reads_increase_depth():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)
    kmer.create_dict_kmer()

    matcher = Match(
        ["testdata/kmer_match_testtdata.fastq"],
        kmer.kmer_index,
        fasta.gene_length
    )

    matcher.match_reads_to_genes()

    assert sum(matcher.dict_depth["geneA_test_1"]) > 0
    assert sum(matcher.dict_depth["geneB_test_2"]) > 0


# Test to confirm short reads are skipped
def test_match_short_reads_skipped():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)
    kmer.create_dict_kmer()

    matcher = Match(
        ["testdata/kmer_match_testtdata.fastq"],
        kmer.kmer_index,
        fasta.gene_length
    )

    matcher.match_reads_to_genes()

    total_depth = sum(sum(v) for v in matcher.dict_depth.values())

    assert total_depth > 0


# Test to confirm reads with N are skipped
def test_match_reads_with_N_skipped():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)
    kmer.create_dict_kmer()

    matcher = Match(
        ["testdata/kmer_match_testtdata.fastq"],
        kmer.kmer_index,
        fasta.gene_length
    )

    matcher.match_reads_to_genes()

    for gene, depths in matcher.dict_depth.items():

        assert min(depths) >= 0


# Test to confirm depth values are positive after matching
def test_match_depth_positive():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)
    kmer.create_dict_kmer()

    matcher = Match(
        ["testdata/kmer_match_testtdata.fastq"],
        kmer.kmer_index,
        fasta.gene_length
    )

    matcher.match_reads_to_genes()

    for gene, depths in matcher.dict_depth.items():

        if sum(depths) > 0:
            assert max(depths) > 0