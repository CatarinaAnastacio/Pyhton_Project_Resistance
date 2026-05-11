"""
Unit tests for the Analysis class.

The Analysis class filters and deduplicates detected resistance genes:
  - calculate_coverage_depth: keep genes that meet minimum coverage and depth
  - join_genes_by_similarity: collapse near-identical alleles
"""

import pytest

from src.analysis import Analysis


##### calculate_coverage_depth #####

# Test to confirm fully covered genes pass the filter
def test_analysis_full_coverage_pass():

    gene_length = {"geneX": 60}
    dict_depth = {"geneX": [50] * 60}

    analysis = Analysis(dict_depth, gene_length, {}, {"geneX": 0})
    result = analysis.calculate_coverage_depth()

    assert "geneX" in result
    assert result["geneX"]["coverage"] == 100.0
    assert result["geneX"]["depth"] == 50.0


# Test to confirm low depth genes are filtered out
def test_analysis_low_depth_removed():

    gene_length = {"geneX": 60}
    dict_depth = {"geneX": [1] * 60}

    analysis = Analysis(dict_depth, gene_length, {}, {"geneX": 0})
    result = analysis.calculate_coverage_depth()

    assert "geneX" not in result


# Test to confirm low coverage genes are filtered out
def test_analysis_low_coverage_removed():

    gene_length = {"geneX": 60}
    dict_depth = {"geneX": [50] * 30 + [0] * 30}

    analysis = Analysis(dict_depth, gene_length, {}, {"geneX": 0})
    result = analysis.calculate_coverage_depth()

    assert "geneX" not in result


# Test to confirm short genes are skipped
def test_analysis_short_gene_skipped():

    gene_length = {"tinyGene": 30}
    dict_depth = {"tinyGene": [50] * 30}

    analysis = Analysis(dict_depth, gene_length, {}, {"tinyGene": 0})
    result = analysis.calculate_coverage_depth()

    assert "tinyGene" not in result


# Test to confirm min depth is calculated correctly
def test_analysis_min_depth():

    gene_length = {"geneX": 60}
    core = [50] * 12 + [11] + [50] * 11
    dict_depth = {"geneX": [0] * 18 + core + [0] * 18}

    analysis = Analysis(dict_depth, gene_length, {}, {"geneX": 0})
    result = analysis.calculate_coverage_depth()

    assert result["geneX"]["min_depth"] == 11




##### join_genes_by_similarity #####

# Test to confirm distinct genes are both kept
def test_analysis_distinct_genes_kept():

    gene_length = {"geneA": 60, "geneB": 60}
    resistance_kmer = {
        "AAAAAAAAAAAAAAAAAAA": [("geneA", 0, "+")],
        "CCCCCCCCCCCCCCCCCCC": [("geneB", 0, "+")],
    }
    filtered = {
        "geneA": {"coverage": 100, "depth": 50, "min_depth": 40},
        "geneB": {"coverage": 100, "depth": 50, "min_depth": 40},
    }

    analysis = Analysis({}, gene_length, resistance_kmer,
                        {"geneA": 10, "geneB": 10})
    final, removed = analysis.join_genes_by_similarity(filtered)

    assert removed == set()
    assert "geneA" in final
    assert "geneB" in final


# Test to confirm similar genes keep the winner with best coverage
def test_analysis_best_coverage_kept():

    gene_length = {"geneA": 60, "geneB": 60}
    resistance_kmer = {}

    for i in range(5):
        resistance_kmer[f"AAAAAAAAAAAAAAAAA{i:02d}"] = [
            ("geneA", i, "+"),
            ("geneB", i, "+"),
        ]
    filtered = {
        "geneA": {"coverage": 100, "depth": 50, "min_depth": 40},
        "geneB": {"coverage": 96,  "depth": 50, "min_depth": 40},
    }

    analysis = Analysis({}, gene_length, resistance_kmer,
                        {"geneA": 10, "geneB": 10})
    final, removed = analysis.join_genes_by_similarity(filtered)

    assert "geneA" in final
    assert "geneB" in removed


# Test to confirm similar genes keep the winner with best depth
def test_analysis_best_depth_kept():

    gene_length = {"geneA": 60, "geneB": 60}
    resistance_kmer = {}

    for i in range(5):
        resistance_kmer[f"AAAAAAAAAAAAAAAAA{i:02d}"] = [
            ("geneA", i, "+"),
            ("geneB", i, "+"),
        ]
    filtered = {
        "geneA": {"coverage": 100, "depth": 80, "min_depth": 40},
        "geneB": {"coverage": 100, "depth": 50, "min_depth": 40},
    }

    analysis = Analysis({}, gene_length, resistance_kmer,
                        {"geneA": 10, "geneB": 10})
    final, removed = analysis.join_genes_by_similarity(filtered)

    assert "geneA" in final
    assert "geneB" in removed


# Test to confirm ties are broken using unique hits
def test_analysis_unique_hits_break_tie():

    gene_length = {"geneA": 60, "geneB": 60}
    resistance_kmer = {}
    
    for i in range(5):
        resistance_kmer[f"AAAAAAAAAAAAAAAAA{i:02d}"] = [
            ("geneA", i, "+"),
            ("geneB", i, "+"),
        ]
    filtered = {
        "geneA": {"coverage": 100, "depth": 50, "min_depth": 40},
        "geneB": {"coverage": 100, "depth": 50, "min_depth": 40},
    }

    analysis = Analysis({}, gene_length, resistance_kmer,
                        {"geneA": 100, "geneB": 10})
    final, removed = analysis.join_genes_by_similarity(filtered)

    assert "geneA" in final
    assert "geneB" in removed
