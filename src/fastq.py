"""
FASTQ reader
Commented version in the full code
"""

import gzip


class Fastq:
    """Reads sequencing reads from a (possibly gzipped) FASTQ file.
    """

    def __init__(self, filename):
        self.filename = filename

    def fastqread(self):
        """Generator that yields just the sequence string of every read."""
        if self.filename.endswith('.gz'):
            infile = gzip.open(self.filename, 'rt')
        else:
            infile = open(self.filename, 'r')

        while True:
            header = infile.readline()
            if not header:
                break

            seq = infile.readline().rstrip()
            infile.readline()  
            infile.readline()   

            if len(seq) == 0:
                continue

            yield seq

        infile.close()
