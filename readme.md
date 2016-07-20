# **SlideSort**

C++ implementation of SlideSort

# Summary
This library provides the implementation of the SlideSort and a related clustering tool [1].

* [1]   Shimizu K, Tsuda K, SlideSort: All pairs similarity search for short reads. Bioinformatics, 2011: 27 (4): 464-470.

# Installation

  git clone https://github.com/iskana/SlideSort.git
  cd SlideSort
  /bin/sh mkall_v1.sh
  /bin/sh mkall_v2.sh

# Quick start with a sample file

  cd SlideSort/SS_v2
  ./slidesort_v2 -d 3 -i test.fasta -mt 4
  or
  cat test.fasta| ./slidesort_v2 -d 3 -I -mt 4
  
  cd SlideSort/mst_v2
  ./ssmst_v2 -d 3 -i test.fasta -mt 4
  or
  cat test.fasta| ./ssmst_v2 -d 3 -I -mt 4 -Z

# Copyright Notice

Copyright (C) 2015, Kana Shimizu
All rights reserved.

# License

See the readme included in this repository.

# Contact Information

* Kana Shimizu (shimizu.kana@waseda.jp)

# History

2016/July/20; initial upload

