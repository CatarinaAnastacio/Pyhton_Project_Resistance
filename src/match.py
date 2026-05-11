"""
Match: align reads to reference genes through shared k-mers.
Commented version in the full code
"""

from src.fastq import Fastq


class Match:
    """Maps each read to the gene it most likely came from.
    """

    def __init__(self, fastq_files, resistance_kmer, gene_lengths):
        self.fastq_files = fastq_files
        self.resistance_kmer = resistance_kmer
        self.gene_lengths = gene_lengths

        self.dict_depth = {gene: [0] * length for gene, length in gene_lengths.items()}
        self.unique_match_counts = {gene: 0 for gene in gene_lengths}

    def match_reads_to_genes(self):
        """Main matching loop -> Iterates over every read in every FASTQ"""
        for fastq_file in self.fastq_files:
            fastq = Fastq(fastq_file)

            for read in fastq.fastqread():
                read_len = len(read)

                if read_len < 19:
                    continue

                possible_hit = False
                for pos in range(0, read_len - 18, 20):
                    kmer = read[pos:pos + 19]
                    if 'N' in kmer:
                        continue
                    if kmer in self.resistance_kmer:
                        possible_hit = True
                        break
                if not possible_hit:
                    continue

                reads_found_gene = {}

                for read_pos in range(read_len - 18):
                    kmer_read = read[read_pos:read_pos + 19]
                    if 'N' in kmer_read or kmer_read not in self.resistance_kmer:
                        continue

                    genes_matched = self.resistance_kmer[kmer_read]

                    first_gene_matched = genes_matched[0][0]
                    unique = True
                    for entry in genes_matched:
                        if entry[0] != first_gene_matched:
                            unique = False
                            break

                    for gene, gene_pos, strand in genes_matched:
                        alignment_shift = gene_pos - read_pos
                        id = (gene, alignment_shift, strand)

                        if id not in reads_found_gene:
                            reads_found_gene[id] = [0, 0]

                        reads_found_gene[id][0] += 1
                        if unique:
                            reads_found_gene[id][1] += 1

                if not reads_found_gene:
                    continue

                best_hits = 0
                for info in reads_found_gene.values():
                    if info[0] > best_hits:
                        best_hits = info[0]

                if best_hits < 3:
                    continue

                for read_pos in range(read_len - 18):
                    kmer_read = read[read_pos:read_pos + 19]
                    if 'N' in kmer_read:
                        continue

                    genes_matched = self.resistance_kmer.get(kmer_read)
                    if genes_matched is None:
                        continue

                    for gene, gene_pos, strand in genes_matched:
                        alignment_shift = gene_pos - read_pos
                        id = (gene, alignment_shift, strand)

                        if id in reads_found_gene and \
                           reads_found_gene[id][0] == best_hits:
                            gene_len = self.gene_lengths[gene]

                            if strand == '+':
                                start = alignment_shift + read_pos
                            else:
                                start = gene_len - 19 - (alignment_shift + read_pos)

                            if 0 <= start <= gene_len - 19:
                                for p in range(start, start + 19):
                                    self.dict_depth[gene][p] += 1

                for (gene, shift, strand), (total_hits, unique_hits) in \
                        reads_found_gene.items():
                    if total_hits == best_hits:
                        self.unique_match_counts[gene] += unique_hits
