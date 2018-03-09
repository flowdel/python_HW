import pandas as pd
import argparse
from Bio import SeqIO

class Kmer_info:


    count = 0
    sequence = ''
    pos = ''

    def __init__(self, kmer_name):
        self.sequence = kmer_name

    def counter(self):
        self.count += 1

    def position(self, index):
        if self.pos != '':
            self.pos += ', '
            self.pos += str(index+1)
        else:
            self.pos += str(index+1)


def find_kmers(file, kmer_size):
    
    sequence = SeqIO.read(file, 'fasta')
    sequence = sequence.seq[0:100]


    kmer_dict = {}

    seq_lng = len(sequence)
    
    for index in range(seq_lng-kmer_size+1):
        current_kmer = sequence[index:(index+kmer_size)]

        if current_kmer in kmer_dict:
            kmer_dict[current_kmer].counter()
            kmer_dict[current_kmer].position(index)
        else:
            kmer_dict[current_kmer] = Kmer_info(current_kmer)
            kmer_dict[current_kmer].counter()
            kmer_dict[current_kmer].position(index)


    kmer_df = pd.DataFrame([])
    for current_kmer in kmer_dict.keys():
        kmer_df = kmer_df.append(pd.DataFrame({'kmer': str(kmer_dict[current_kmer].sequence),
                                               'number': kmer_dict[current_kmer].count,
                                               'positions': str(kmer_dict[current_kmer].pos)}, index=[0]), ignore_index=True)


    sorted_kmer_df = kmer_df.sort_values(['number'], ascending=False)
    print(sorted_kmer_df.iloc[[0]])

if __name__ == "__main__":
       
    parser = argparse.ArgumentParser(description='Find kmers')
    
    parser.add_argument('-f', '--file', help='sequence', type=str, required=True)
    parser.add_argument('-s', '--kmer_size', help='kmer size', type=int, default=3)
    
    args = parser.parse_args()
    file = args.file
    kmer_size = args.kmer_size
    
    find_kmers(file, kmer_size)


