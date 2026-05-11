"""
Tests for the Fasta class.

The reference test FASTA 'gene_kmer_dict.fsa' contains two genes of
length 60 with very repetitive sequences
"""

import pytest

from src.fasta import Fasta


# Test to confirm FASTA reading returns the correct number of genes
def test_fasta_number_of_genes():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")

    genes = list(fasta.reading_dnaFasta())

    assert len(genes) == 2


# Test to confirm headers are read correctly
def test_fasta_headers():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")

    headers = [header for header, _ in fasta.reading_dnaFasta()]

    assert "geneA_test_1" in headers
    assert "geneB_test_2" in headers


# Test to confirm sequence lengths are correct
def test_fasta_sequence_lengths():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")

    for _, sequence in fasta.reading_dnaFasta():

        assert len(sequence) == 60


# Test to confirm gene lengths are stored correctly
def test_fasta_gene_lengths():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")

    fasta.read_fasta_dict()

    assert fasta.gene_length["geneA_test_1"] == 60
    assert fasta.gene_length["geneB_test_2"] == 60


# Test to confirm gene sequences are stored correctly
def test_fasta_gene_sequences():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")

    fasta.read_fasta_dict()

    expected = "ACGT" * 15

    assert fasta.gene_sequences["geneA_test_1"] == expected


# Test to confirm both dictionaries have the same keys
def test_fasta_dictionary_keys():

    fasta = Fasta("testdata/gene_kmer_dict.fsa")

    fasta.read_fasta_dict()

    assert set(fasta.gene_length.keys()) == set(fasta.gene_sequences.keys())


# Test to confirm nonexistent files raise an exception
def test_fasta_nonexistent_file():

    with pytest.raises(Exception):

        fasta = Fasta("testdata/notexisting.fsa")

        list(fasta.reading_dnaFasta())