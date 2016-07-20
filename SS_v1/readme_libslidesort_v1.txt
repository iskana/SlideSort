/**************************************************/
                   LIBSLIDESORT (ver.1)
    Copyright (c) 2010,   Kana Shimizu
              All rights reserved.
/**************************************************/

LIBSLIDESORT is a library of SLIDESORT.
Similar pairs are obtained from a user defined call-back function.

1. How to use library

Just link libslidesort when you compile your source.

1. 1 Sample source code.

The sample source code (test_libslidesort.cpp) calculates degree of input sequences.


compile:
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./
gcc -o test_libslidesort test_libslidesort.cpp -L. -lslidesort

execution:
./test_libslidesort -i test.fasta -d 3 -l


2. How to obtain options

Options can be obtained from command line arguments in the same manner as binary version of SLIDESORT.

*** I RECOMMEND TO USE -l OPTION which enables the program to avoid outputting pairs. This option is very useful for a large-scale dataset. ****


3. Dependency of shared libraries for Linux OS
ss-intel-64:
        libm.so.6
        libstdc++.so.6
        libgcc_s.so.1
        libc.so.6
        libdl.so.2
        /lib64/ld-linux-x86-64.so.2
ss-gcc-64:
        libm.so.6
        libc.so.6
        /lib64/ld-linux-x86-64.so.2

4. Contact

shimizu.kana@waseda.jp


5. Reference

Kana Shimizu and Koji Tsuda "SlideSort:All pairs similarity search for short reads", Bioinformatics, 2011: 27(4): 464-470.

6. Lisence

LIBSLIDESORT is freely available both for academic and commercial user.
Lisencing is necessary if LIBSLIDESORT is redistributed for commercial purpose.
Redistributions in any form must reproduce the above copyright notice.
