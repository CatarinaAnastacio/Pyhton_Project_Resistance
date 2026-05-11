"""
FASTA reader for the resistance gene reference database
Commented version in the full code
"""


class Fasta:
    """Reads a FASTA file with the reference resistance genes
    """

    def __init__(self, filename):
        self.filename = filename
        self.gene_length = {}
        self.gene_sequences = {}

    def reading_dnaFasta(self):
        """Generator that yields (header, sequence) for every gene recors
        """
        infile = open(self.filename, 'r')
        header = ''
        sequence = ''

        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if header != '':
                    yield header, sequence
                header = line[1:]   # drop the '>'
                sequence = ''
            else:
                sequence += line

        if sequence != '':
            yield header, sequence

        infile.close()

    def read_fasta_dict(self):
        """Iterate the file once and fill the two dictionaries"""
        for header, seq in self.reading_dnaFasta():
            name = header.strip().split(' ')[0]
            self.gene_length[name] = len(seq)
            self.gene_sequences[name] = seq
