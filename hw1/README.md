# SEQUENCE ALIGNMENT

This program will allow you to make pairwise alignment of sequences that are in the same fasta file.

## Description

Program was created on Python, and it uses global alignment algorithm. 

### Required arguments

* -f - file with seqiences in fasta format
* -b - you can use Blosum matrix or not
* -m - match value
* -ms - mismatch value
* -g - gap value

### Running

You can use test dataset to check program's work.

Open terminal and run:
```
sequence_alignment.py -f test_data.fasta -b 
```
or

```
sequence_alignment.py -f test_data.fasta -m 5 -ms -3 -g -10
```

## Addition

* [BLOSUM](https://en.wikipedia.org/wiki/BLOSUM) - You can read about BLOSUM here


## Author

* **Adel Gazizova** - *BI student* - [flowdel](https://github.com/flowdel)
