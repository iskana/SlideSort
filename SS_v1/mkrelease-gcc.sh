#!/bin/sh

rm -rf ss_v1-gcc-64
rm -f ss_v1-gcc-64.tar.gz

chmod 755 mdcommon.h.perl

cp common.h common.h.org
make -f Makefile.gcc clean

./mdcommon.h.perl common.h.org  VCPP SSE4 CALL_BACK_FUNC INTERNAL_GAP_OPEN CALL_BACK_FUNC_BY_VIRTUAL_FUNC> common.h
make -f Makefile.gcc

make -f Makefile.gcc clean

./mdcommon.h.perl common.h.org  VCPP SSE4 INTERNAL_GAP_OPEN > common.h
make -f Makefile.gcc slib

mkdir ss_v1-gcc-64

cp slidesort_v1 libslidesort_v1.so common.h param.h mscls.h test.fasta readme_v1.txt test_libslidesort_v1.cpp readme_libslidesort_v1.txt ./ss_v1-gcc-64
echo g++ -o test_libslidesort_v1 test_libslidesort_v1.cpp -L. -lslidesort_v1 > ./ss_v1-gcc-64/make.sh

make -f Makefile.gcc clean
make -f Makefile.gcc lib

tar cvf ss_v1-gcc-64.tar ss_v1-gcc-64
gzip ss_v1-gcc-64.tar
cp libslidesort_v1.a ./ss_v1-gcc-64/

cp common.h.org common.h
