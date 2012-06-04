#!/bin/sh
## code by H.Canary
# $1 is the install path which overwrites the default.
# $2 is the build directory path.

if test "$1" = "-h" ; then
	echo "Useage"
	echo "  ${0} [INSTALL_PREFIX] [BUILD_DIR]"
	exit 1
fi

if test "$1" ; then
	INSTALL_PREFIX="$1"
	mkdir -p "$INSTALL_PREFIX"
	INSTALL_PREFIX=`cd "$INSTALL_PREFIX" ; pwd`
else
	INSTALL_PREFIX="${HOME}/local"
	mkdir -p "$INSTALL_PREFIX"
fi

if test "$2" ; then
	BUILD_DIR="$2"
	mkdir -p "$BUILD_DIR"
else
	BUILD_DIR="build_emulator"
	mkdir -p "$BUILD_DIR"
fi

SRC_DIR=`dirname "$0"`
SRC_DIR=`cd "$SRC_DIR" ; pwd`

echo "SRC_DIR=\"${SRC_DIR}\""
echo "BUILD_DIR=\"${BUILD_DIR}\""
echo "INSTALL_PREFIX=\"${INSTALL_PREFIX}\""

cd "$BUILD_DIR"

cmake \
	-D "CMAKE_INSTALL_PREFIX:PATH=${INSTALL_PREFIX}" \
	"$SRC_DIR" && \
	make && \
	make install

{	echo '#!/bin/sh'
	echo "LD_LIBRARY_PATH=\"${INSTALL_PREFIX}/lib\""
	echo "export LD_LIBRARY_PATH"
	echo "exec \"${INSTALL_PREFIX}/bin/emulator\" \"\$@\""
} > "${INSTALL_PREFIX}/bin/run_emulator"
chmod +x "${INSTALL_PREFIX}/bin/run_emulator"

{	echo '#!/bin/sh'
	echo "LD_LIBRARY_PATH=\"${INSTALL_PREFIX}/lib\""
	echo "export LD_LIBRARY_PATH"
	echo "exec \"${INSTALL_PREFIX}/bin/estimator\" \"\$@\""
} > "${INSTALL_PREFIX}/bin/run_estimator"
chmod +x "${INSTALL_PREFIX}/bin/run_estimator"
