#!/bin/sh

if [ -e /opt/intel/Compiler/11.1/056/bin/intel64/ ]; then
flag=intel
else
flag=gcc
fi

cd ./SS_v1

cp comlineoption.cpp comlineoption.org
cat comlineoption.cpp |grep -v "mt number of logical processors used for calculation" > tmp.cpp
mv tmp.cpp comlineoption.cpp

/bin/sh ./mkrelease-gcc.sh

if [ $flag = "intel" ] ; then
/bin/sh ./mkrelease-intel.sh
fi

cp comlineoption.org comlineoption.cpp

mv ss_v1-gcc-64.tar.gz ss_v1-gcc-64-LL.tar.gz

if [ $flag = "intel" ] ; then
mv ss_v1-intel-64.tar.gz ss_v1-intel-64-LL.tar.gz 
mv ss_v1-intel-SSE4-64.tar.gz  ss_v1-intel-SSE4-64-LL.tar.gz
fi

cd ../mst_v1
/bin/sh ./mkrelease.sh $flag

mv ssmst_v1-gcc-64.tar.gz ssmst_v1-gcc-64-LL.tar.gz

if [ $flag = "intel" ] ; then
mv ssmst_v1-intel-64.tar.gz ssmst_v1-intel-64-LL.tar.gz
mv ssmst_v1-intel-SSE4-64.tar.gz ssmst_v1-intel-SSE4-64-LL.tar.gz
fi
