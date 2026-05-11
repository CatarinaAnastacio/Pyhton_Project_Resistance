"""
Tests for the Kmer class.

Both genes in 'gene_kmer_dict.fsa' are 60 bases long, so:
    - 60 - 18 = 42 forward k-mers per gene
    - 42 reverse-complement k-mers per gene
    - 2 genes  -> 168 total entries
"""

import pytest

from src.fasta import Fasta
from src.kmer import Kmer
from src.dna import DNA


# Test to confirm that the kmer index is created
def test_kmer_index_not_empty():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)

    kmer.create_dict_kmer()

    assert len(kmer.kmer_index) > 0


# Test to confirm the total number of kmer entries
def test_kmer_total_entries():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)

    kmer.create_dict_kmer()

    total_entries = sum(len(v) for v in kmer.kmer_index.values())

    assert total_entries == 168


# Test to confirm that each entry has the correct format
def test_kmer_entry_shape():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)

    kmer.create_dict_kmer()

    for entries in kmer.kmer_index.values():

        for gene, position, strand in entries:

            assert isinstance(gene, str)
            assert isinstance(position, int)
            assert strand in ('+', '-')


# Test to confirm the first forward kmer of geneA is present
def test_first_forward_kmer_present():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)

    kmer.create_dict_kmer()

    first_kmer = "ACGTACGTACGTACGTACG"

    assert first_kmer in kmer.kmer_index

    entries = kmer.kmer_index[first_kmer]

    assert ("geneA_test_1", 0, '+') in entries


# Test to confirm the reverse complement kmer is present
def test_reverse_complement_kmer_present():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")
    fasta.read_fasta_dict()

    kmer = Kmer(fasta)

    kmer.create_dict_kmer()

    reverse_kmer = DNA().reverse_com_strand("ACGTACGTACGTACGTACG")

    assert reverse_kmer in kmer.kmer_index

    entries = kmer.kmer_index[reverse_kmer]

    strands = [strand for _, _, strand in entries]

    assert '-' in strands
