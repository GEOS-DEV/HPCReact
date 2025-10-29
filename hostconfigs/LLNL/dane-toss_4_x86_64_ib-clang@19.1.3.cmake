#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: clang@=19.1.3
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-19.1.3-magic/bin/clang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-19.1.3-magic/bin/clang++" CACHE PATH "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "")

#--------------------------------------------------------------------------------
# CMake Standard
#--------------------------------------------------------------------------------

set(BLT_CXX_STD "c++17" CACHE STRING "")


#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_MKL ON CACHE BOOL "")

set(MKL_INCLUDE_DIRS "/usr/tce/packages/mkl/mkl-2022.1.0/include" CACHE PATH "")

set(MKL_LIBRARIES /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_intel_lp64.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_gnu_thread.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_core.so
                  /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-13.3.1/llvm-19.1.3-gy2lu5xbi4csr2k47emlajzfs5mlsd4g/bin/../lib/x86_64-unknown-linux-gnu/libomp.so
                  /lib64/libpthread.so
                  /lib64/libm.so
                  /lib64/libdl.so CACHE STRING "")