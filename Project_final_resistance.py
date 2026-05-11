'''
Antibiotic Resistance Gene Detection in Metagenomic NGS Data
============================================================

Goal: Having a FASTA file of resistance genes and one or more gzipped FASTQ
files containing metagenomic sequencing reads, the goal is to identify which resistance
genes are present in the sample.

Method (high level):
    1. Build a 19-mer index from the resistance gene database (both strands).
    2. Go through every FASTQ read and look for 19-mers that hit the index.
    3. For reads with several k-mers that align on the same gene at consistent
       positions (which means a real alignment), increase the depth of those positions.
    4. Filter genes by coverage (>= 95%) and average depth (>= 10).
    5. Resolution of near-identical homologs.

IMPORTANT: The full print of statistical features and plots of final genes and their homologs is in GitHub: https://github.com/CatarinaAnastacio/Pyhton_Project_Resistance
'''



import gzip


#################
##### FASTA #####
#################

class Fasta:
    '''
    Reads a FASTA file with the reference resistance genes.

    Computing of two dictionaries with gene name as:
        gene_length[name] -> length of the gene in bases
        gene_sequences[name] -> the DNA sequence as a string
    '''
    def __init__(self, filename):
        self.filename = filename
        self.gene_length = {}
        self.gene_sequences = {}

    def reading_dnaFasta(self):
        '''
        Generator that yields (header, sequence) for every gene in the file.

        Using yield avoids keeping the whole file in memory, which is
        good if the database is large.
        '''
        infile = open(self.filename, 'r')
        header = ''
        sequence = ''

        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                # New record, finalize the last one  
                if header != '':
                    yield header, sequence
                header = line[1:]  # don´t keep the '>' character
                sequence = ''
            else:
                # sequence line
                sequence += line
        if sequence != '':
            yield header, sequence

        infile.close()

    def read_fasta_dict(self):
        '''
        Iterate the file once and fill the two dictionaries
        Returns gene_length and gene_sequences with the gene as key
        '''
        for header, seq in self.reading_dnaFasta():
            name = header.strip().split(' ')[0]  #  keep only the first part of header (the gene name)
            self.gene_length[name] = len(seq)
            self.gene_sequences[name] = seq



###############
##### DNA #####
###############

class DNA:
    '''
    DNA function: reverse complement and 19-mer extraction
    '''
    def reverse_com_strand(self, sequence):
        '''
        Returns the reverse complement of the input DNA sequence

        Reads can come from both DNA strands, so we need both
        orientations of every reference gene
        '''
        comp = {'A':'T','T':'A','G':'C','C':'G'}
        rev = []
        # Go through the sequence from the end and complement each base
        for base in sequence[::-1]:
            if base in comp:
                rev.append(comp[base])
            else:
                # Anything not recognised becomes an 'N'
                rev.append('N')
        return ''.join(rev)

    def nineteen_pos(self, DNA_seq):
        '''
        Generator that yields (position, kmer) for every 19-mer in DNA_seq
        K-mers that contain 'N' are skipped because
        they would give unreliable matches
        '''
        for i in range(len(DNA_seq) - 18):
            kmer = DNA_seq[i:i+19]
            if 'N' not in kmer:
                yield i, kmer

################
##### KMER #####
################
class Kmer:
    '''
    Creates a hash table aligning every reference 19-mer to where it happens

    Structure:
        kmer_index[kmer] = [(gene_name, position, strand), ...]

    Both strands of every gene are indexed.
    '''
    def __init__(self, fasta):
        self.fasta = fasta
        self.kmer_index = {}
        self.dna = DNA()

    def create_dict_kmer(self):
        '''
        Computes the kmers for all gene sequences, both forward and reverse.
        Each k-mer maps to a list of (gene_name, position, strand) entries
        '''
        for header, sequence in self.fasta.reading_dnaFasta():
            gene_name = header.strip().split(' ')[0]

            # Forwar strand
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


#################
##### FASTQ #####
#################

class Fastq:
    '''
    Reads sequencing reads from a gzipped FASTQ file

    A FASTQ record has 4 lines: @header, sequence, +, quality
    We only need the sequence line, so we drop the rest
    '''
    def __init__(self, filename):
        self.filename = filename

    def fastqread(self):
        '''
        Generator that yields the sequence of each read in the file
        '''
        if self.filename.endswith('.gz'):
            infile = gzip.open(self.filename, 'rt')
        else:
            infile = open(self.filename, 'r')

        while True:
            header = infile.readline()
            if not header:
                break
            seq = infile.readline().rstrip()
            infile.readline()   # the '+' line
            infile.readline()   # the quality line

            if len(seq) == 0:
                break

            yield seq

        infile.close()


####################
##### MATCHING #####
####################

class Match:
    '''
    Maps each read to the gene from where that read came from

    For every read:
      1. Doa pre-filter (check a few k-mers and skip reads with no match)
      2. For all read k-mers, look them up in the reference index.
      3. Group hits by (gene, alignment_shift, strand). Hits sharing the
         same shift come from a good local alignment. Isolated hits are not good.
      4. Keep only the best alignment for that read (with most matching k-mers).
      5. Update depth per position for the gene at every covered position
    '''

    def __init__(self, fastq_files, resistance_kmer, gene_lengths):
        self.fastq_files = fastq_files
        self.resistance_kmer = resistance_kmer
        self.gene_lengths = gene_lengths

        # depth[gene][pos] = how many reads cover position 'pos' of the gene
        self.dict_depth = {gene: [0] * length for gene, length in gene_lengths.items()}

        # count read k-mers that only match one gene
        self.unique_match_counts = {gene: 0 for gene in gene_lengths}

    def match_reads_to_genes(self):
        '''
        Matching loop. Iterates over every read in every FASTQ
        Map sequencing reads to reference genes using k-mer matches
        '''
        for fastq_file in self.fastq_files:
            fastq = Fastq(fastq_file)

            for read in fastq.fastqread():
                read_len = len(read)

                # Reads shorter than k can't have any 19-mers
                if read_len < 19:
                    continue
                
                # Pre-filter: only one in 20 k-mers is sampled
                # Most reads in a sample don't come from resistance genes, so this saves time
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


                # Go through every k-mer of the read and group hits by (gene, alignment_shift, strand)
                # If the read came from this gene, the shift
                # would map every read base to the corresponding gene base.

                # reads_found_gene[id] = [total_hits, unique_hits]
                reads_found_gene = {}

                for read_pos in range(read_len - 18):
                    kmer_read = read[read_pos:read_pos + 19]
                    if 'N' in kmer_read or kmer_read not in self.resistance_kmer:
                        continue

                    genes_matched = self.resistance_kmer[kmer_read]

                    # If the kmer is specific to a single gene, we count them as 'unique' to that gene
                    # This will count as evidence after to filter
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


                # Choose the best alignment for the current read, which is the 
                # (gene, shift, strand) with the most matching k-mers. 
                # Reject reads with less than 3 k-mer hits, because
                # 1 or 2 hits are probably sequencing errors 
                best_hits = 0
                for info in reads_found_gene.values():
                    if info[0] > best_hits:
                        best_hits = info[0]

                if best_hits < 3:
                    continue
                

                # Update depth: walk the read again and increment every gene position
                # covered by a matching k-mer. Do it only for alignments tied with best_hits count,
                # because one read can equally support several homologous genes.
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

                        if id in reads_found_gene and reads_found_gene[id][0] == best_hits:
                            gene_len = self.gene_lengths[gene]

                            # Turn k-mer position in read to the position in the gene
                            # For the reverse strand we flip the position
                            # back onto the forward strand of the gene
                            if strand == '+':
                                start = alignment_shift + read_pos
                            else:
                                start = gene_len - 19 - (alignment_shift + read_pos)

                            if 0 <= start <= gene_len - 19:
                                for p in range(start, start + 19):
                                    self.dict_depth[gene][p] += 1

                # Count of unique k-mer hits per gene
                # this is then used to break ties between similar genes.
                for (gene, shift, strand), (total_hits, unique_hits) in reads_found_gene.items():
                    if total_hits == best_hits:
                        self.unique_match_counts[gene] += unique_hits

##########################
##### ANALYSIS CLASS #####
##########################

class Analysis:
    '''
    Computes coverage / depth and resolves homologous genes.

    Only 'core region' is used because there is a natural drop in
    depth at the gene ends (due to incomplete k-mer coverage) 
    '''
    def __init__(self, dict_depth, gene_length, resistance_kmer, dict_unique_hits):
        self.dict_depth = dict_depth
        self.gene_length = gene_length
        self.resistance_kmer = resistance_kmer
        self.dict_unique_hits = dict_unique_hits

    def calculate_coverage_depth(self, min_cov=95, min_depth=10):
        '''Filter genes by minimum core coverage and average core depth'''
        filtered_genes = {}

        for gene, pos_counts in self.dict_depth.items():
            gene_len = self.gene_length[gene]
            
            # Skip the first and last 18 bases as they are not be fully covered by
            # full 19-mers, so they would lower coverage
            core_start = 18
            core_end = gene_len - 18

            if core_end <= core_start: 
                # gene is too short for the analysis
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
        '''
        Collapse near-identical genes

        Two genes are considered the same 'family' if they share at least
        the 'threshold' of the smaller gene's k-mer set. In a 'family' we
        keep the gene with the strongest evidence (coverage, then depth,
        then number of unique k-mer hits if there is a tie).
        '''

        # Build a kmer set per gene for intersection
        gene_kmers = {gene: set() for gene in self.gene_length}
        for kmer, entries in self.resistance_kmer.items():
            for gene, position, strand in entries:
                gene_kmers[gene].add(kmer)

        gene_list = list(filtered_genes.keys())
        removed = set()

        # Comparing all (after filtering is not that heavy)
        for i in range(len(gene_list)):
            for j in range(i + 1, len(gene_list)):
                gene_a = gene_list[i]
                gene_b = gene_list[j]

                if gene_a in removed or gene_b in removed:
                    continue

                kmer_a = gene_kmers[gene_a]
                kmer_b = gene_kmers[gene_b]

                shared = len(kmer_a.intersection(kmer_b))   # computing shared kmers

                if shared == 0:
                    continue

                smaller = kmer_a if len(kmer_a) < len(kmer_b) else kmer_b

                similarity = shared / len(smaller)

                if similarity >= threshold:
                    # they are homologs; now we device which one wins
                    stat_a = filtered_genes[gene_a]
                    stat_b = filtered_genes[gene_b]

                    if stat_a['coverage'] != stat_b['coverage']:
                        better, worse = (gene_a, gene_b) if stat_a['coverage'] > stat_b['coverage'] else (gene_b, gene_a)
                    elif stat_a['depth'] != stat_b['depth']:
                        better, worse = (gene_a, gene_b) if stat_a['depth'] > stat_b['depth'] else (gene_b, gene_a)
                    else:
                        # In case of a tie, check number of unique k-mer matches
                        better, worse = (gene_a, gene_b) if self.dict_unique_hits[gene_a] >= self.dict_unique_hits[gene_b] else (gene_b, gene_a)

                    removed.add(worse)

        final = {}
        for gene in filtered_genes:
            if gene not in removed:
                final[gene] = filtered_genes[gene]


        return final, removed
    

# ----------------------------
# MAIN
# ----------------------------
if __name__ == '__main__':

    # Input files (edit to match other data)
    fasta_file = 'resistance_genes.fsa.txt'
    fastq_files = ['Unknown3_raw_reads_1.txt.gz', 'Unknown3_raw_reads_2.txt.gz']

    # Load reference gene database
    fasta = Fasta(fasta_file)
    fasta.read_fasta_dict()

    # Build k-mer index
    kmer = Kmer(fasta)
    kmer.create_dict_kmer()

    # Match every read against the index
    matcher = Match(fastq_files, kmer.kmer_index, fasta.gene_length)
    matcher.match_reads_to_genes()

    # Compute coverage and depth; Collapse homologs
    analysis = Analysis(matcher.dict_depth, fasta.gene_length,
                        kmer.kmer_index, matcher.unique_match_counts)

    filtered_genes = analysis.calculate_coverage_depth()
    final_genes, removed = analysis.join_genes_by_similarity(filtered_genes)

    # Sort and Print
    # Sort: highest coverage first, then highest avg depth, then highest min depth
    sorted_genes = sorted(final_genes.items(),
        key=lambda x: (-x[1]['coverage'], -x[1]['depth'], -x[1]['min_depth']))
    
    print('Resistance genes present in the metagenomic samples:')
    for gene, stats in sorted_genes: 
        print(gene, "| coverage:", f"{stats['coverage']:.1f}",
                    "| depth:", f"{stats['depth']:.3f}",
                    "| min_depth:", round(stats['min_depth'])
                )