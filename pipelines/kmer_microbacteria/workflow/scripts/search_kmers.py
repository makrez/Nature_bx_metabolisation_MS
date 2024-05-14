# Input file is a fasta file with the kmers
from Bio import SeqIO
import fire
import os
import pandas as pd

# annotated_fasta = '/data/projects/p533_Plant_root_bacteria_genomes/' +\
#     'genome_assembly_pipeline/results/4_prokka/LMB2/LMB2.ffn'
#
# kmer_file = '/data/projects/p533_Plant_root_bacteria_genomes/kmer_analysis/' +\
#     'results/03_kmerscores/score9.fasta'

def extract_matches(annotated_fasta, kmer_file, output_file):
    matches = list()
    for fasta in SeqIO.parse(annotated_fasta, 'fasta'):
        samplename = os.path.basename(annotated_fasta).strip('.ffn')
        for kmer in SeqIO.parse(kmer_file, 'fasta'):
            kmers_in_seq = list()
            revcomp = kmer.seq.reverse_complement()

            kmer_in_fasta = fasta.seq.count(kmer.seq)
            revcomp_in_fasta = fasta.seq.count(revcomp)

            if kmer_in_fasta + revcomp_in_fasta > 0:
                matches.append(
                        [samplename, fasta.description, str(kmer.seq), kmer.name]
                    )

    m = pd.DataFrame(matches,
                     columns = ['Sample', 'Description', 'kmer', 'kmer_name'])

    with open(output_file, 'w') as csv:
        m.to_csv(csv, index=False)

if __name__ == '__main__':
    fire.Fire(extract_matches)
