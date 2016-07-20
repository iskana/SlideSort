/**************************************************/
                   LIBSLIDESORT (ver.2)
    Copyright (c) 2010,   Kana Shimizu
              All rights reserved.
/**************************************************/

LIBSLIDESORT is a library of SLIDESORT.
Similar pairs are obtained from a user defined call-back function.

1. How to use library

Just link libslidesort when you compile your source.

1. 1 Sample source code.

The sample source code (test_libslidesort_v2.cpp) calculates degree of input sequences.


How to compile:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./
g++ -L. -lslidesort_v2 -o test_libslidesort_v2 test_libslidesort_v2.cpp -fopenmp

An execution example:
./test_libslidesort -i test.fasta -d 3 -l -mt 3


2. How to obtain options

Options can be obtained from command line arguments in the same manner as binary version of SLIDESORT.

*** I RECOMMEND TO USE -l OPTION which enables the program to avoid outputting pairs. This option is very useful for a large-scale dataset. ****


3. Contact

shimizu.kana@waseda.jp


4. Reference

Kana Shimizu and Koji Tsuda "SlideSort:All pairs similarity search for short reads", Bioinformatics, 2011: 27(4): 464-470.

5. Lisence

LIBSLIDESORT is freely available both for academic and commercial user.
Lisencing is necessary if LIBSLIDESORT is redistributed for commercial purpose.
Redistributions in any form must reproduce the above copyright notice.

