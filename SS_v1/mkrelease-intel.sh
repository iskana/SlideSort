#!/bin/sh

rm -rf ss_v1-intel-SSE4-64
rm -f  ss_v1-intel-SSE4-64.tar.gz 

rm -rf ss_v1-intel-64
rm -f  ss_v1-intel-64.tar.gz 

export PATH=/opt/intel/Compiler/11.1/056/bin/intel64/:$PATH
chmod 755 mdcommon.h.perl

mkdir ss_v1-src
cp ckExecution.sh ./ss_v1-src
cp Makefile.gcc ./ss_v1-src
cp Makefile.intel ./ss_v1-src
cp mdcommon.h.perl ./ss_v1-src
cp mkrelease-gcc.sh ./ss_v1-src
cp mkrelease-intel.sh ./ss_v1-src
cp readme_v1.txt ./ss_v1-src
cp readme_libslidesort_v1.txt ./ss_v1-src
cp test.fasta ./ss_v1-src
cp test_libslidesort_v1.cpp ./ss_v1-src
cp *.cpp ./ss_v1-src
cp *.h ./ss_v1-src
tar cvf ss_v1-src.tar ./ss_v1-src
gzip ss_v1-src.tar

make -f Makefile.intel clean

cp common.h common.h.org

./mdcommon.h.perl common.h.org  CALL_BACK_FUNC INTERNAL_GAP_OPEN  CALL_BACK_FUNC_BY_VIRTUAL_FUNC> common.h
make -f Makefile.intel

make -f Makefile.intel clean

./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN > common.h
make -f Makefile.intel slib

mkdir ss_v1-intel-SSE4-64

cp slidesort_v1 libslidesort_v1.so common.h param.h mscls.h test.fasta readme_v1.txt test_libslidesort_v1.cpp readme_libslidesort_v1.txt ./ss_v1-intel-SSE4-64
echo icpc -o test_libslidesort_v1 test_libslidesort_v1.cpp -L. -lslidesort_v1 > ./ss_v1-intel-SSE4-64/make.sh

make -f Makefile.intel clean
make -f Makefile.intel lib

tar cvf ss_v1-intel-SSE4-64.tar ss_v1-intel-SSE4-64
gzip ss_v1-intel-SSE4-64.tar
cp libslidesort_v1.a ./ss_v1-intel-SSE4-64/

make -f Makefile.intel clean

./mdcommon.h.perl common.h.org  CALL_BACK_FUNC INTERNAL_GAP_OPEN SSE4  CALL_BACK_FUNC_BY_VIRTUAL_FUNC> common.h
make -f Makefile.intel

make -f Makefile.intel clean

./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN SSE4 > common.h
make -f Makefile.intel slib

mkdir ss_v1-intel-64

cp slidesort_v1 libslidesort_v1.so common.h param.h mscls.h test.fasta readme_v1.txt test_libslidesort_v1.cpp readme_libslidesort_v1.txt ./ss_v1-intel-64
echo icpc -o test_libslidesort_v1 test_libslidesort_v1.cpp -L. -lslidesort_v1 > ./ss_v1-intel-64/make.sh

make -f Makefile.intel clean
make -f Makefile.intel lib

tar cvf ss_v1-intel-64.tar ss_v1-intel-64
gzip ss_v1-intel-64.tar
cp libslidesort_v1.a ./ss_v1-intel-64/

cp common.h.org common.h
