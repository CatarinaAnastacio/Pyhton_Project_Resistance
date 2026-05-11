"""
Analysis: turns depth per position into a final gene list with coverage and avg depth.
Commented version in the full code
"""


class Analysis:
    """Computes coverage and depth, and resolves homologous genes.
    """

    def __init__(self, dict_depth, gene_length, resistance_kmer, dict_unique_hits):
        self.dict_depth = dict_depth
        self.gene_length = gene_length
        self.resistance_kmer = resistance_kmer
        self.dict_unique_hits = dict_unique_hits

    def calculate_coverage_depth(self, min_cov=95, min_depth=10):
        """Filter genes by minimum core coverage and average core depth"""
        filtered_genes = {}

        for gene, pos_counts in self.dict_depth.items():
            gene_len = self.gene_length[gene]

            core_start = 18
            core_end = gene_len - 18

            if core_end <= core_start:
                continue

            core = pos_counts[core_start:core_end]
            n = len(core)

            covered = sum(1 for pos_depth in core if pos_depth > 0)
            coverage = covered / n * 100
            avg_depth = sum(core) / n

            if coverage >= min_cov and avg_depth >= min_depth:
                filtered_genes[gene] = {
                    'coverage': coverage,
                    'depth': avg_depth,
                    'min_depth': min(core) if core else 0
                }

        return filtered_genes

    def join_genes_by_similarity(self, filtered_genes, threshold=0.5):
        """Collapse near-identical variants
        """

        gene_kmers = {gene: set() for gene in self.gene_length}
        for kmer, entries in self.resistance_kmer.items():
            for gene, position, strand in entries:
                gene_kmers[gene].add(kmer)

        gene_list = list(filtered_genes.keys())
        removed = set()

        for i in range(len(gene_list)):
            for j in range(i + 1, len(gene_list)):
                gene_a = gene_list[i]
                gene_b = gene_list[j]

                if gene_a in removed or gene_b in removed:
                    continue

                kmer_a = gene_kmers[gene_a]
                kmer_b = gene_kmers[gene_b]

                shared = len(kmer_a.intersection(kmer_b))
                if shared == 0:
                    continue

                smaller = kmer_a if len(kmer_a) < len(kmer_b) else kmer_b
                similarity = shared / len(smaller)

                if similarity >= threshold:
                    stat_a = filtered_genes[gene_a]
                    stat_b = filtered_genes[gene_b]

                    if stat_a['coverage'] != stat_b['coverage']:
                        better, worse = (gene_a, gene_b) if stat_a['coverage'] > stat_b['coverage'] \
                            else (gene_b, gene_a)
                    elif stat_a['depth'] != stat_b['depth']:
                        better, worse = (gene_a, gene_b) if stat_a['depth'] > stat_b['depth'] \
                            else (gene_b, gene_a)
                    else:
                        better, worse = (gene_a, gene_b) \
                            if self.dict_unique_hits[gene_a] >= self.dict_unique_hits[gene_b] \
                            else (gene_b, gene_a)

                    removed.add(worse)

        final = {gene: stats for gene, stats in filtered_genes.items() if gene not in removed}
        return final, removed
