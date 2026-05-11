# Resistance to Antibiotics — Project 22118

## Project structure

```
project_root/
 ── Project_final_resistance.py         # full version of the code 
 ── src/                                # version in classes
     ── __init__.py
     ── fasta.py                 # FASTA reader
     ── dna.py                   # reverse complement + 19-mer extraction
     ── kmer.py                  # builds the k-mer hash index
     ── fastq.py                 # gzipped FASTQ reader
     ── match.py                 # read-to-gene matching algorithm
     ── analysis.py              # coverage/depth filtering + homolog resolution
 ── tests/                       # 42 pytest unit tests
     ── test_fasta.py
     ── test_dna.py
     ── test_kmer.py
     ── test_fastq.py
     ── test_match.py
     ── test_analysis.py
 ── testdata/                    # test data used by the tests
     ── gene_kmer_dict.fsa
     ── kmer_match_testtdata.fastq.gz
     ── empty_reads.fastq.gz
     ── no_match_reads.fastq.gz 
```

## Running the pipeline

Edit the file paths in the end of `Project_final_resistance.py` to run other data:

```bash
python Project_final_resistance.py
```

The program prints all genes that pass the 95 % coverage and 10 depth thresholds,
sorted by coverage, average depth, then minimum depth.

## Running the tests

From the project path:

```bash
pip install pytest
pytest test/ -v
```

All 42 tests should pass.
