/**************************************************/
                   SLIDESORT (ver.2)
    Copyright (c) 2010-2024,   Kana Shimizu
              All rights reserved.
/**************************************************/

SLIDESORT is a fast and exact method that can find all similar pairs
from a string pool in terms of edit distance.

The input to SLIDESORT is a set of sequences of almost equal length
and a threshold of edit distance.
From the input sequences, SLIDESORT exactly finds all pairs within
the input edit-distance threshold.
SLIDESORT also accepts a gap opening cost and a gap extension cost.

SLIDESORT can find similar pairs
  - within edit-distance d, 
  - from sequences whose length range from minimum length L to L+d 
  - with arbitrary gap open cost G and gap extension cost g.

1. Input format

Multi-fasta format and FASTQ format are acceptable.
(see http://en.wikipedia.org/wiki/FASTA_format)
SLIDESORT accepts both DNA sequence and protein sequence.

ex)
> seq 1
ATGCTAGCTGATACATCTAGCTCGTACGTACGTCAGTCGTAGT
CTGACTGACTAGCTAGCTAGCATCGTACGTACGTCGTAGCTAC
> seq 2
ATGCTAGCTGATACATCTAGCTCGTACGTACGTCAGTCGTAGT
CTGACTGACTAGCTAGCTAGCATCGTACGTACGTCGTAGCTAC


2. Execution

2.1 typical usage

./slidesort -d <distance>  -i <input_fasta_file> -o <output_file>

2.2 Finding pairs from a part of the input sequences

./slidesort -p -d <edit-distance-threshold> -fst_head <starting_line> -fst_size <size_of_sequences> -i <input_file> -o <output_file>

2.3 Comparing two datasets

./slidesort -m -d <edit-distance-threshold> -fst_head <starting_line_of_dataset1> -fst_size <size_of_sequences_of_dataset1>
-snd_head <starting_line_of_dataset2> -snd_size <size_of_sequences_of_dataset2> -i <input_file> -o <output_file>

2.4 Sequences including unknown characters.
In default setting,  SlideSort excludes sequences with unknown characters such as 'N', 'X'.
To include sequences with unknown characters, use "-u" option.

3 Options
slidesort -d <distance>  [options]

help: 
  -h show this message.

Important options:
  -d  distance threshold
  -i  input filename  (default=input.txt)
  -o  output filename  (default=stdout)
  -t  distance type  E: edit-distance  H: hamming-distance (default=E)
  -v  search with both original seq and reverse complement seq.
  -mt number of logical processors used for calculation.
Advanced options:
  -a  output alignment
  -c  type of input string.
      DNA: DNA seq, PROTEIN: protein seq, INT: integer seq (default=DNA)
  -g  gap extention cost (default=1, must be positive value. it is better to use larger value to avoid slow-down of the search.)
  -G  gap open cost (default=0, must be positive value)
  -k  size of sorting key
  -l  do not output results. (for library version)
  -m  use cross search mode. (find pairs between two datasets. first set: from fst_head to fst_head+fst_size. second set snd_head to  
snd_head+snd_size)
  -p  use partial mode. (find pairs from fst_head to fst_head+fst_size)
  -u  do not exclude sequences with unknown character. ex) n, N, Z, etc...
  -V  output a same pair twice if dist(A,B)<=d and dist(A,B')<=d. B' is reverse complement of B. (use with -v)\n\


4. Library (LibSlideSort)
See readme_libslidesort.txt.

5. Contact

shimizu.kana@waseda.jp


6. Reference

Kana Shimizu and Koji Tsuda "SlideSort:All pairs similarity search for short reads", Bioinformatics, 2011: 27(4): 464-470.

7. Lisence

SLIDESORT is freely available both for academic and commercial users.
Redistributions in binary form must reproduce the above copyright notice.
