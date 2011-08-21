#!/bin/sh
# $1 is the install path which overwritest the one given
# in CMakeLists.txt
cd build

if test "$1" 
then
		# strangely this doesn't seem to work, at least not on os-x
		echo "using install prefix: $1"
		cmake .. -DCMAKE_INSTALL_PREFIX=${1}
else
		cmake ..
fi
make && make install

cd ..