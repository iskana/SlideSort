#!/bin/sh

export PATH=/opt/intel/Compiler/11.1/056/bin/intel64/:$PATH

rm -rf ss_v2-gcc-64
rm -r ss_v2-gcc-64.tar.gz
rm -rf ss_v2-intel-SSE4-64
rm -r ss_v2-intel-SSE4-64.tar.gz
rm -rf ss_v2-intel-SSE4-64
rm -r ss_v2-intel-SSE4-64.tar.gz

cp ../SS_v1/mdcommon.h.perl ./
cp ../SS_v1/*.cpp ./
cp ../SS_v1/*.h ./

echo "#include \"mscls.h\"" >tmp.h
cat parallelslidesort.h |grep -v "mscls.h" >>tmp.h
mv tmp.h parallelslidesort.h

chmod 755 mdcommon.h.perl
cp common.h common.h.org
make -f Makefile.gcc clean

echo "COMPILING MSLIDESORT BY GCC"

echo "MAKE STAND ALONE PROGRAM"
./mdcommon.h.perl common.h.org  VCPP SSE4 INTERNAL_GAP_OPEN > common.h
make -f Makefile.gcc

make -f Makefile.gcc clean

echo "MAKE SHARED LIBRARY"
make -f Makefile.gcc slib

mkdir ss_v2-gcc-64
cp slidesort_v2 libslidesort_v2.so common.h param.h parallelslidesort.h mscls.h test.fasta readme_v2.txt test_libslidesort_v2.cpp readme_libslidesort_v2.txt ./ss_v2-gcc-64
echo g++ -o test_libslidesort_v2 test_libslidesort_v2.cpp -L. -lslidesort_v2 -fopenmp > ./ss_v2-gcc-64/make.sh

make -f Makefile.gcc clean

echo "MAKE STATIC LIBRARY"
make -f Makefile.gcc lib

tar cvf ss_v2-gcc-64.tar ss_v2-gcc-64
gzip ss_v2-gcc-64.tar
cp libslidesort_v2.a ./ss_v2-gcc-64/

cp common.h.org common.h

if [ $1 = "intel" ] ; then
echo "COMPILING MSLIDESORT BY INTEL COMPILER"

./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN > common.h
make -f Makefile.intel

make -f Makefile.intel clean

echo "MAKE SHARED LIBRARY"
make -f Makefile.intel slib

mkdir ss_v2-intel-SSE4-64
cp slidesort_v2 libslidesort_v2.so common.h param.h parallelslidesort.h mscls.h test.fasta readme_v2.txt test_libslidesort_v2.cpp readme_libslidesort_v2.txt ./ss_v2-intel-SSE4-64
echo icpc -o test_libslidesort_v2 test_libslidesort_v2.cpp -L. -lslidesort_v2 -openmp > ./ss_v2-intel-SSE4-64/make.sh

make -f Makefile.intel clean

echo "MAKE STATIC LIBRARY"
make -f Makefile.intel lib

tar cvf ss_v2-intel-SSE4-64.tar ss_v2-intel-SSE4-64
gzip ss_v2-intel-SSE4-64.tar
cp libslidesort_v2.a ./ss_v2-intel-SSE4-64/

cp common.h.org common.h
fi

