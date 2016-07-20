/**************************************************/
                    SSMST
    Copyright (c) 2010-2012,   Kana Shimizu
              All rights reserved.
/**************************************************/

SSMST is a fast tool that can find minimum spanning trees from a string pool.
Distance between two strings are evaluated either by edit-distance or hamming-distance.

The input to SSMST is a set of sequences of almost equal length
and a threshold of edit-distance/hamming-distance.

1. Input format

Multi-fasta format and fastq are acceptable.
(see http://en.wikipedia.org/wiki/FASTA_format)
SSMST accepts both DNA sequence and protein sequence.

ex)
> seq 1
ATGCTAGCTGATACATCTAGCTCGTACGTACGTCAGTCGTAGT
CTGACTGACTAGCTAGCTAGCATCGTACGTACGTCGTAGCTAC
> seq 2
ATGCTAGCTGATACATCTAGCTCGTACGTACGTCAGTCGTAGT
CTGACTGACTAGCTAGCTAGCATCGTACGTACGTCGTAGCTAC


2. Execution

2.1 typical usage

./ssmst -d <distance>  -i <input_fasta_file> -o <output_file>

2.2 Options
ssmst -d <distance>  [options]

help:
  -h show this message.

Important options:
  -d  distance threshold
  -t  distance type  E: edit-distance  H: hamming-distance (default=E)
  -i  input filename  (default=input.txt)
  -o  output filename  (default=tree.txt)
  -v  search with both original seq and reverse complement seq.
Advanced options:
  -T  type of tree format. (S:SIF, O:SlideSort original)
  -c  type of input string.
      DNA: DNA seq, PROTEIN: protein seq, INT: integer seq (default=DNA)
  -u  do not exclude sequence with unknown character. ex) n, N, Z, etc...
  -g  gap extention cost (default=1, same value as mismatch cost)
  -G  gap open cost (default=0, must be positive value)
  -k  size of sorting key
  -p  use partial mode. (find pairs from fst_head to fst_head+fst_size)
  -m  use cross search mode. (find pairs between two datasets. first set: from fst_head to fst_head+fst_size. second set snd_head to snd_head+snd_size)
  -V  output a same pair twice if dist(A,B)<=d and dist(A,B')<=d. B' is reverse complement of B. (use with -v)

3. Output
MSTs are described in following two format.

(1) SSMST original format
 - Trees are separated by "--".
 - A left most node of each line is a parent node, and other nodes are it's children.
 - Distance is shown in ().

For examples of the following trees, 

tree_1:
A
+-----+--+
B     C  D
+--+  |  +--+
E  F  G  H  I

dist(A,B)=1
dist(A,C)=2
dist(A,D)=1
dist(D,H)=1
dist(D,I)=1
dist(C,G)=2
dist(B,E)=1
dist(B,F)=1

id:tree_1--
A: B(1), C(2), D(1),
D: H(1), I(1),
C: G(2),
B: E(1), F(1),

(2) SIF (simple interaction format) 
seq1 - seq2	distance
seq2 - seq3	distance
...

3. Contact

shimizu.kana@waseda.jp


4. Reference

Kana Shimizu and Koji Tsuda "SlideSort:All pairs similarity search for short reads", Bioinformatics, 2011: 27(4): 464-470.

5. Lisence
SSMST is freely available both for academic and commercial user.
Lisencing is necessary if SSMST is redistributed for commercial purpose.
Redistributions in any form must reproduce the above copyright notice.



