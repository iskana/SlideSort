#!/bin/sh

rm -rf ssmst_v2-gcc-64
rm -rf ssmst_v2-intel-64
rm -rf ssmst_v2-intel-SSE4-64
rm -rf ssmst_v2-gcc-64.tar.gz
rm -rf ssmst_v2-intel-64.tar.gz
rm -rf ssmst_v2-intel-SSE4-64.tar.gz

rm ssmst_v2

export PATH=/opt/intel/Compiler/11.1/056/bin/intel64/:$PATH
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
chmod 755 mdcommon.h.perl

cp ../SS_v2/mdcommon.h.perl ./


cp ../SS_v2/*.h ./
cp ../SS_v2/*.cpp ./
cp ../SS_v1/*.h ./
cp ../SS_v1/*.cpp ./
cp ../mst_v1/olmst.h ./
cp ../mst_v1/mstree.cpp ./
cp ../mst_v1/tree.cpp ./


echo "#include \"parallelslidesort.h\"" >tmp.h
cat olmst.h |grep -v "parallelslidesort.h" >>tmp.h
mv tmp.h olmst.h

cp common.h common.h.org

rm ssmst_v2
make -f  Makefile.gcc clean
./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN VCPP SSE4 SINGLE_CORE > common.h
make -f  Makefile.gcc 
mkdir ssmst_v2-gcc-64
cp ssmst_v2 ssmst_v2-gcc-64
cp readme-ssmst_v2.txt ssmst_v2-gcc-64
make -f  Makefile.gcc clean

if [ $1 = "intel" ] ; then
./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN SINGLE_CORE > common.h
make -f  Makefile.intel 
mkdir ssmst_v2-intel-SSE4-64
cp ssmst_v2 ssmst_v2-intel-SSE4-64
cp readme-ssmst_v2.txt ssmst_v2-intel-SSE4-64

cp common.h.org common.h
fi

tar cvf ssmst_v2-gcc-64.tar ssmst_v2-gcc-64
gzip ssmst_v2-gcc-64.tar

if [ $1 = "intel" ] ; then
tar cvf ssmst_v2-intel-SSE4-64.tar ssmst_v2-intel-SSE4-64
gzip ssmst_v2-intel-SSE4-64.tar
fi