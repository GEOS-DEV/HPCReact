cmake_minimum_required( VERSION 3.23.1 )

# Set version number
set( HPCREACT_VERSION_MAJOR 0 )
set( HPCREACT_VERSION_MINOR 1 )
set( HPCREACT_VERSION_PATCHLEVEL 0 )

# check if this is build as a submodule or a separate project
get_directory_property( parent_dir PARENT_DIRECTORY )
if(parent_dir)
    set( is_submodule ON )
else()
    set( is_submodule OFF )
endif()

if( NOT is_submodule )
    message( "not a submodule")
    project( hpcReact LANGUAGES CXX C )

    set( BLT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/cmake/blt CACHE PATH "" )
    set( BLT_CXX_STD "c++17" CACHE STRING "Version of C++ standard" FORCE )
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)

    set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "" )

    option( HPCREACT_ENABLE_UNIT_TESTS "Builds tests" ON )

    option( ENABLE_CUDA "Build with CUDA" OFF )
    option( ENABLE_HIP "Build with HIP" OFF )

endif()

include( ${BLT_SOURCE_DIR}/SetupBLT.cmake )

include( cmake/CMakeBasics.cmake )
include( cmake/Macros.cmake )
#include( cmake/Config.cmake )


add_subdirectory( src )

# if( HPCREACT_ENABLE_DOCS )
#     add_subdirectory( docs )
# endif()


