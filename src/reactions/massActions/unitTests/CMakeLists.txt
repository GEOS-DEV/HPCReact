# Specify list of tests
set( testSourceFiles
     testMassActions.cpp
     )

set( dependencyList hpcReact gtest )

# Add gtest C++ based tests
foreach(test ${testSourceFiles})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )
    blt_add_test( NAME ${test_name}
                   COMMAND ${test_name} )
endforeach()