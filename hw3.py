from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import argparse

class Kmer_spectrum:
       
    def __init__(self):
        self.file = ''
        self.sequence = ''
        self.file = ''
        self.kmer_dict = {}
        self.counter = 0
        self.spectrum_data = []
        self.max_yvalue = 0
    
    
    def kmer_search(self, file, kmer_size, q):
        
        for record in SeqIO.parse(file, 'fastq'):
            ok = True
            qual = record.letter_annotations['phred_quality']
            for i in qual:
                if i < q:
                    ok = False
            if ok:
                read = record.seq
                read_lng = len(record.seq)
                for index in range(read_lng-kmer_size+1):
                    current_kmer = read[index:(index+kmer_size)]
                    if current_kmer in self.kmer_dict:
                        self.kmer_dict[current_kmer] = self.kmer_dict[current_kmer] + 1
                    else:
                        self.kmer_dict[current_kmer] = 1
                        
            else:
                continue

    def spectrum(self):
        
        self.spectrum_data = Counter(self.kmer_dict.values()).most_common()
        lst_y = []
        #print(self.spectrum_data)
        for i in Counter(self.kmer_dict.values()).most_common():
            lst_y.append(i[1])
        self.max_yvalue = int(max(lst_y)/200)
        #print(self.max_yvalue)
        count, frequency = zip(*self.spectrum_data)
        plt.bar(count, frequency, width=1, color='red')
        plt.ylim((0, self.max_yvalue))
        plt.tight_layout()
        plt.show() 
    
    def spectrum_min(self, minimum):
        
        count, frequency = zip(*self.spectrum_data)
        fig, ax = plt.subplots()
        ax.bar(count, frequency, width=1, color='red')
        plt.ylim((0, self.max_yvalue))
        plt.axvline(x=int(minimum))
        plt.tight_layout()
        plt.show() 
    
    def genome_size(self, minimum):
        c_1 = 0
        c_2 = 0
        c_3 = 0
        for i in range(len(self.spectrum_data)):
            if self.spectrum_data[i][0] > int(minimum):
                c_1 += self.spectrum_data[i][0]*self.spectrum_data[i][1]
                c_2 += self.spectrum_data[i][0]
                c_3 += 1
        c_2 = c_2/c_3   
        size = c_1/c_2
        print("Approximate genome size: ", int(size))
        
if __name__ == "__main__":
       
    parser = argparse.ArgumentParser(description='Kmer spectrum')
    
    parser.add_argument('-f', '--file', help='fastq file', type=str, required=True)
    parser.add_argument('-k', '--kmer_size', help='kmer size', type=int, default=15)
    parser.add_argument('-q', '--quality', help='min quality of each nucleotide in fatsq file', type=int, default=20)
    
    args = parser.parse_args()
    file = args.file
    kmer_size = args.kmer_size
    q = args.quality
    
    data = Kmer_spectrum()
    data.kmer_search(file, kmer_size, q)
    #data.kmer_search('/Volumes/Element/Копия_test_kmer.fastq', 15, 20)
    data.spectrum()
    minimum = int(input('Please enter the noise cut-off boundary: '))
    data.spectrum_min(minimum)
    data.genome_size(minimum)