message( "this hostconfig assumes you are using homebrew")
message( "brew install bison cmake gcc git-lfs open-mpi openblas python3 llvm cppcheck lcov")

message( "CMAKE_SYSTEM_PROCESSOR = ${CMAKE_SYSTEM_PROCESSOR}" )
message("CONFIG_NAME = ${CONFIG_NAME}")

set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI OFF CACHE PATH "")
# set(MPI_C_COMPILER "${HOMEBREW_DIR}/bin/mpicc" CACHE PATH "")
# set(MPI_CXX_COMPILER "${HOMEBREW_DIR}/bin/mpicxx" CACHE PATH "")
# set(MPIEXEC "${HOMEBREW_DIR}/bin/mpirun" CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)

set(ENABLE_CUDA "OFF" CACHE PATH "" FORCE)
set(ENABLE_OPENMP "OFF" CACHE PATH "" FORCE)

set(ENABLE_CALIPER "OFF" CACHE PATH "" FORCE )

set(ENABLE_DOXYGEN ON CACHE BOOL "" FORCE)
set( DOXYGEN_EXECUTABLE "${HOMEBREW_DIR}/bin/doxygen" CACHE PATH "" FORCE)

set( CPPCHECK_EXECUTABLE "${HOMEBREW_DIR}/bin/cppcheck" CACHE PATH "" FORCE)
set( CLANGTIDY_EXECUTABLE "${HOMEBREW_DIR}/opt/llvm/bin/clang-tidy" CACHE PATH "" FORCE)

set( ENABLE_COVERAGE OFF CACHE BOOL "" FORCE)
set( GCOV_EXECUTABLE "/usr/bin/gcov" CACHE PATH "" FORCE)
set( LCOV_EXECUTABLE "${HOMEBREW_DIR}/bin/lcov" CACHE PATH "" FORCE)


