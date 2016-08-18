#!/bin/sh

rm -rf ssmst_v1-gcc-64
rm -rf ssmst_v1-intel-64
rm -rf ssmst_v1-intel-SSE4-64
rm -rf ssmst_v1-gcc-64.tar.gz
rm -rf ssmst_v1-intel-64.tar.gz
rm -rf ssmst_v1-intel-SSE4-64.tar.gz

rm -rf ssmst_v1

export PATH=/opt/intel/Compiler/11.1/056/bin/intel64/:$PATH
export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH
chmod 755 mdcommon.h.perl

cp ../SS_v1/common.h ./
cp ../SS_v1/param.h ./
cp ../SS_v1/mscls.h ./
cp ../SS_v1/mdcommon.h.perl ./

cp mstree.cpp mstree.org
cat mstree.cpp |grep -v "mt number of logical processors used for calculation" > tmp.cpp
mv tmp.cpp mstree.cpp

cp olmst.h olmst.h.org
echo "#include \"mscls.h\"" >tmp.h
cat olmst.h |grep -v "mscls.h">>tmp.h
mv tmp.h olmst.h

cp common.h common.h.org

if [ $1 = "intel" ] ; then
cp ../SS_v1/ss_v1-intel-SSE4-64/libslidesort_v1.a ./
./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN > common.h
make -f Makefile.intel
mkdir ssmst_v1-intel-SSE4-64
cp ssmst_v1 ssmst_v1-intel-SSE4-64
cp readme-ssmst_v1.txt ssmst_v1-intel-SSE4-64

rm ssmst_v1
make -f Makefile.intel clean
cp ../SS_v1/ss_v1-intel-64/libslidesort_v1.a ./
./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN SSE4 > common.h
make -f Makefile.intel
mkdir ssmst_v1-intel-64
cp ssmst_v1 ssmst_v1-intel-64
cp readme-ssmst_v1.txt ssmst_v1-intel-64
fi

rm ssmst_v1
make -f Makefile.intel clean
cp ../SS_v1/ss_v1-gcc-64/libslidesort_v1.a ./
./mdcommon.h.perl common.h.org  INTERNAL_GAP_OPEN VCPP SSE4 > common.h
make -f Makefile.gcc 
mkdir ssmst_v1-gcc-64
cp ssmst_v1 ssmst_v1-gcc-64
cp readme-ssmst_v1.txt ssmst_v1-gcc-64

cp common.h.org common.h
cp mstree.org mstree.cpp

tar cvf ssmst_v1-gcc-64.tar ssmst_v1-gcc-64
gzip ssmst_v1-gcc-64.tar

if [ $1 = "intel" ] ; then
tar cvf ssmst_v1-intel-64.tar ssmst_v1-intel-64
gzip ssmst_v1-intel-64.tar
tar cvf ssmst_v1-intel-SSE4-64.tar ssmst_v1-intel-SSE4-64
gzip ssmst_v1-intel-SSE4-64.tar
fi

mv olmst.h.org olmst.h

