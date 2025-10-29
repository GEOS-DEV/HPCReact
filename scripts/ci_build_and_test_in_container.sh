#!/bin/bash
env


if [ -f /etc/os-release ]; then
  . /etc/os-release
  if [[ "$ID" == "ubuntu" ]]; then
      PACKAGE_MANAGER="apt"
  elif [[ "$ID" == "rocky" ]]; then
      PACKAGE_MANAGER="yum"
  else
      echo "Unsupported OS: $ID"
      exit 1
  fi
else
  echo "/etc/os-release not found. Unable to determine OS."
  exit 1
fi


if [ -f /usr/lib64/libblas.so.3 ]; then
  BLAS_LIB=/usr/lib64/libblas.so.3
  LAPACK_LIB=/usr/lib64/liblapack.so.3
elif [ -f /usr/lib/x86_64-linux-gnu/libblas.so ]; then
  BLAS_LIB=/usr/lib/x86_64-linux-gnu/libblas.so
  LAPACK_LIB=/usr/lib/x86_64-linux-gnu/liblapack.so
else
  echo "BLAS/LAPACK not found"; 

  echo "Using package manager: $PACKAGE_MANAGER"

  if [ "$PACKAGE_MANAGER" == "apt" ]; then
      apt update && apt-get install -y libblas-dev liblapack-dev
  elif [ "$PACKAGE_MANAGER" == "yum" ]; then
      yum update && yum install -y blas lapack 
  fi

fi








# The or_die function run the passed command line and
# exits the program in case of non zero error code
function or_die () {
    "$@"
    local status=$?
    echo status = $status
    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}

# Working in the root of the cloned repository
or_die cd $(dirname $0)/..

if [[ -z "${HOST_CONFIG}" ]]; then
  echo "Environment variable \"HOST_CONFIG\" is undefined."
  exit 1
fi

if [[ -z "${CMAKE_BUILD_TYPE}" ]]; then
  echo "Environment variable \"CMAKE_BUILD_TYPE\" is undefined."
  exit 1
fi

if [[ "$*" == *--code-coverage* ]]; then
  ENABLE_COVERAGE=ON
else
  ENABLE_COVERAGE=OFF
fi

HPCREACT_BUILD_DIR=/tmp/build
HPCREACT_INSTALL_DIR=/tmp/install
or_die python3 scripts/config-build.py \
               -hc ${HOST_CONFIG} \
               -bt ${CMAKE_BUILD_TYPE} \
               -bp ${HPCREACT_BUILD_DIR} \
               -ip ${HPCREACT_INSTALL_DIR} \
               -DBLAS_LIBRARIES=${BLAS_LIB} \
               -DLAPACK_LIBRARIES=${LAPACK_LIB} \
               -DENABLE_COVERAGE:BOOL=${ENABLE_COVERAGE}

or_die cd ${HPCREACT_BUILD_DIR}

# Code style check
if [[ "$*" == *--test-code-style* ]]; then
  or_die ctest --output-on-failure -R "testUncrustifyCheck"
  exit 0
fi

# Documentation check
if [[ "$*" == *--test-doxygen* ]]; then

  if [ "$PACKAGE_MANAGER" == "apt" ]; then
      apt update && apt-get install -y texlive-full
  elif [ "$PACKAGE_MANAGER" == "yum" ]; then
      yum update && yum install -y texlive 
  fi

  or_die ctest --output-on-failure -R "testDoxygenCheck"
  exit 0
fi

# code checks
if [[ "$*" == *--code-checks* ]]; then
  or_die ctest --output-on-failure -R "testCppCheck|testClangTidy"
  exit 0
fi



if [[ "$*" == *--code-coverage* ]]; then
  or_die make -j ${NPROC} VERBOSE=1
  or_die make hpcReact_coverage VERBOSE=1
  cp -r ${HPCREACT_BUILD_DIR}/hpcReact_coverage.info.cleaned /tmp/hpcReact/hpcReact_coverage.info.cleaned
fi


if [[ "$*" == *--build-exe* ]]; then
  or_die make -j ${NPROC} VERBOSE=1

  if [[ "$*" != *--disable-unit-tests* ]]; then
    or_die ctest --output-on-failure -E "testUncrustifyCheck|testDoxygenCheck|testCppCheck|testClangTidy"
  fi
fi




exit 0
