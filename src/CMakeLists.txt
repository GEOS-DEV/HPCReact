
set( hpcReact_headers
     reactions/ReactionsBase.hpp
     )

set( hpcReact_sources
     )

blt_add_library( NAME             hpcReact
#                 SOURCES          ${hpcReact_sources}
                 HEADERS          ${hpcReact_headers}
                 DEPENDS_ON       ${hpcReact_dependencies}
                )

# target_include_directories( hpcReact
#                             INTERFACE
#                             $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
#                             $<INSTALL_INTERFACE:include> )

# install( FILES ${hpcReact_headers} 
#          DESTINATION include )


# install( TARGETS hpcReact
#          EXPORT hpcReact-targets
#          ARCHIVE DESTINATION lib
#          LIBRARY DESTINATION lib
#          RUNTIME DESTINATION lib )

# install( EXPORT hpcReact-targets
#          DESTINATION lib/cmake/hpcReact )

# hpcReact_add_code_checks( PREFIX hpcReact
#                        EXCLUDES "blt/*" )

add_subdirectory( reactions/unitTests )