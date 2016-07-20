#!/bin/sh

if [ -e /opt/intel/Compiler/11.1/056/bin/intel64/ ]; then
flag=intel
else
flag=gcc
fi

cd ./SS_v2
/bin/sh ./mkrelease.sh $flag

mv ss_v2-gcc-64.tar.gz ss_v2-gcc-64-LL.tar.gz

if [ $flag = "intel" ] ; then
mv ss_v2-intel-SSE4-64.tar.gz  ss_v2-intel-SSE4-64-LL.tar.gz
fi

cd ../mst_v2
/bin/sh ./mkrelease.sh $flag

mv ssmst_v2-gcc-64.tar.gz ssmst_v2-gcc-64-LL.tar.gz

if [ $flag = "intel" ] ; then
mv ssmst_v2-intel-SSE4-64.tar.gz ssmst_v2-intel-SSE4-64-LL.tar.gz
fi
