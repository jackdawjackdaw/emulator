#!/bin/sh
## code by H.Canary
# $1 is the install path which overwrites the default.
# $2 is the build directory path.

if test "$1" = "-h" ; then
	echo "Useage"
	echo "  ${0} [INSTALL_PREFIX] [BUILD_DIR]"
	exit 1
fi

choosedir() {
        D=${1:-${2}}
        mkdir -p "$D"
        echo `cd "$D" ; pwd`
}

INSTALL_PREFIX=`choosedir "$1" "${HOME}/local"`
BUILD_DIR=`choosedir "$2" "build_test"`

SRC_DIR=`dirname "$0"`
SRC_DIR=`cd "$SRC_DIR" ; pwd`

echo "SRC_DIR=\"${SRC_DIR}\""
echo "BUILD_DIR=\"${BUILD_DIR}\""
echo "INSTALL_PREFIX=\"${INSTALL_PREFIX}\""

cd "$BUILD_DIR"

mkscpt() {
        LIB_PATH="${1}/lib"
        EXE="${1}/bin/$2"
        SCPT="${1}/bin/run_${2}"
        {       echo '#!/bin/sh'
                echo "LD_LIBRARY_PATH=\"${LIB_PATH}\""
                echo "export LD_LIBRARY_PATH"
                echo "exec \"${EXE}\" \"\$@\""
        } > "$SCPT"
        chmod +x "$SCPT"
        echo "created \"${SCPT}\"."
}

cmake \
	-D "CMAKE_INSTALL_PREFIX:PATH=${INSTALL_PREFIX}" \
	"$SRC_DIR" && \
	make && \
	make install #&& \
	# mkscpt "$INSTALL_PREFIX" estimator && \
	# mkscpt "$INSTALL_PREFIX" emulator && \
	# mkscpt "$INSTALL_PREFIX" interactive_emulator
