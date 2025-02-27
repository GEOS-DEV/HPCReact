blt_add_doxygen_target( hpcReact_doxygen )

add_test( NAME testDoxygenCheck
          COMMAND bash -c "${mkdir_cmd} ${DOXYGEN_EXECUTABLE} Doxyfile 2> >(tee doxygen.err) && exit $(cat doxygen.err | wc -l)"
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)