export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/.
export PATH=/opt/intel/Compiler/11.1/056/bin/intel64/:$PATH


cd ./ss-gcc-64/
./slidesort
/bin/sh make.sh
./test_libslidesort

cd ../

cd ./ss-intel-64/
./slidesort
/bin/sh make.sh
./test_libslidesort

cd ../

cd ./ss-intel-SSE4-64/
./slidesort
/bin/sh make.sh
./test_libslidesort

