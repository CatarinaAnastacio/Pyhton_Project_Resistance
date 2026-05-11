"""
DNA function: reverse complement and 19-mer extraction.
Commented version in the full code
"""


class DNA:

    def reverse_com_strand(self, sequence):
        """Return the reverse complement of a DNA sequence.
        """
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        rev = []
        for base in sequence[::-1]:
            if base in comp:
                rev.append(comp[base])
            else:
                rev.append('N')
        return ''.join(rev)

    def nineteen_pos(self, DNA_seq):
        """Generator yielding (position, kmer) for every k-mer (length 19) in DNA_seq.
        """
        for i in range(len(DNA_seq) - 18):
            kmer = DNA_seq[i:i + 19]
            if 'N' not in kmer:
                yield i, kmer
