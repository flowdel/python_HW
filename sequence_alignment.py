import numpy as np
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo as mat
import matplotlib.pyplot as plt
import argparse


blosum = mat.blosum62

def match(pair, blosum):
    if pair in blosum:
        return blosum[pair]
    else:
        return blosum[(tuple(reversed(pair)))]

def score(seq1, seq2, use_blosum, match_value, mismatch_value, gap_value, blosum=blosum): 
    
    if len(seq2) > len(seq1):
        seq1, seq2 = seq2, seq1

    n = len(seq1) + 1
    m = len(seq2) + 1

    Wn = np.zeros((n, m))
    new_m = np.zeros((n, m))
    
    for i in range(1, n):
        new_m[i][0] = 2
    for j in range(1, m):
        new_m[0][j] = 1
    for i in range(1, n):
        Wn[i][0] = (gap_value) * i
    for j in range(1, m):
        Wn[0][j] = (gap_value) * j 
    
    for i in range(1, n):
        for j in range(1, m):
            pair = (seq1[i-1], seq2[j-1])
            a1 = Wn[i][j-1] + gap_value
            a2 = Wn[i-1][j] + gap_value
            if seq1[i-1] != seq2[j-1]:
                if use_blosum:
                    a3 = Wn[i-1][j-1] + match(pair,blosum)
                else:
                    a3 = Wn[i-1][j-1] + mismatch_value  
            else:
                if use_blosum:
                    a3 = Wn[i-1][j-1] + match(pair,blosum)
                else:
                    a3 = Wn[i-1][j-1] + match_value
            Wn[i][j] = max([a1, a2, a3])
            if (a1 > a2) and (a1 > a3):
                new_m[i][j] = 1  # гориз
            elif (a2 > a3):
                new_m[i][j] = 2  # вертик
            else:
                new_m[i][j] = 3  # диаг
    
    sequence1 = []
    sequence2 = []
    p = []
    n = len(seq1)
    m = len(seq2)
    i = n
    j = m
    while new_m[i][j] != 0:
        if new_m[i][j] == 2:
            sequence2.insert(0, '-')
            sequence1.insert(0, seq1[i-1])
            p.insert(0, ' ')
            i -= 1
        elif new_m[i][j] == 1:
            sequence1.insert(0, '-')
            sequence2.insert(0, seq2[j-1])
            p.insert(0, ' ')
            j -= 1
        else:
            sequence1.insert(0, seq1[i-1])
            sequence2.insert(0, seq2[j-1])
            p.insert(0, '|')
            j -= 1
            i -= 1
        
    score = Wn[n-1][m-1]
    return score

def align(file, blosum=blosum):
    
    with open(file, 'r') as f:
        files = SeqIO.parse(f, 'fasta')
        seqs = [fasta.seq for fasta in files]
        m = np.array([[score(x, y, use_blosum, match_value, mismatch_value, gap_value, blosum) for x in seqs] for y in seqs])
        print('Score matrix:')
        print(m)
        plt.imshow(m, interpolation="nearest")
        plt.colorbar(orientation='vertical')

if __name__ == "__main__":
       
    parser = argparse.ArgumentParser(description='Sequences alignment')
    
    parser.add_argument('-f', '--file', help='file with seqiences', type=str, required=True)
    parser.add_argument('-b', '--blosum_using', help='to use or not to use blosum', action='store_true', default=False)
    parser.add_argument('-m', '--match', help='match_value', type=int, default=1)
    parser.add_argument('-ms', '--mismatch', help='mismatch_value', type=int, default=-1)
    parser.add_argument('-g', '--gap', help='gap_value', type=int, default=-5)
    
    args = parser.parse_args()
    sequences = args.file
    use_blosum = args.blosum_using
    match_value = args.match
    mismatch_value = args.mismatch
    gap_value = args.gap
    
    align(sequences, blosum=blosum)
