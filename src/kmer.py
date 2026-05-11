"""
K-mer index for the resistance gene database
Commented version in the full code
"""

from src.dna import DNA


class Kmer:
    """Builds a hash table mapping every reference 19-mer to where it occurs
    """

    def __init__(self, fasta):
        self.fasta = fasta
        self.kmer_index = {}
        self.dna = DNA()

    def create_dict_kmer(self):
        """Goes through every gene and store all 19-mers from + and - strands"""
        for header, sequence in self.fasta.reading_dnaFasta():
            gene_name = header.strip().split(' ')[0]

            # Forward strand
            for position, kmer in self.dna.nineteen_pos(sequence):
                if kmer not in self.kmer_index:
                    self.kmer_index[kmer] = []
                self.kmer_index[kmer].append((gene_name, position, '+'))

            # Reverse complement strand
            rev = self.dna.reverse_com_strand(sequence)
            for rev_position, kmer in self.dna.nineteen_pos(rev):
                if kmer not in self.kmer_index:
                    self.kmer_index[kmer] = []
                self.kmer_index[kmer].append((gene_name, rev_position, '-'))
