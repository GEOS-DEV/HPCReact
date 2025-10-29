set(CMAKE_ENABLE_EXPORTS ON)

if( CMAKE_BUILD_TYPE MATCHES "Debug" )

else()

endif()

option( ENABLE_COVERAGE "Enable coverage" OFF )
#set( SHIVA_BUILD_OBJ_LIBS OFF CACHE BOOL "" )


# if( CMAKE_CXX_STANDARD IN_LIST "98; 11; 14" )
#     MESSAGE(FATAL_ERROR "Shiva requires at least c++17")
# endif()


blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT "${OpenMP_CXX_FLAGS}")
blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS
                                 GNU   "-Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wcast-align -Wcast-qual"
                                 CLANG "-Wpedantic -pedantic-errors -Wshadow -Wfloat-equal -Wcast-align -Wcast-qual -Wabsolute-value"
                               )

blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS_DEBUG
                                 GNU ""
                                 CLANG "-fstandalone-debug"
                                )

if( ENABLE_COVERAGE )
  blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS_DEBUG
                                   GNU   "-O0 -fno-inline -fno-inline-functions -fno-inline-functions-called-once"
                                  )
endif()

set( CAMP_ENABLE_TESTS OFF CACHE BOOL "")
