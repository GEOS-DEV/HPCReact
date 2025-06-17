
set( hpcReact_headers
     common/macros.hpp
     common/CArrayWrapper.hpp
     reactions/bulkGeneric/EquilibriumReactions.hpp
     reactions/bulkGeneric/EquilibriumReactionsAggregatePrimaryConcentration_impl.hpp
     reactions/bulkGeneric/EquilibriumReactionsReactionExtents_impl.hpp
     reactions/geochemistry/GeochemicalSystems.hpp
     reactions/bulkGeneric/KineticReactions.hpp
     reactions/bulkGeneric/KineticReactions_impl.hpp
     reactions/bulkGeneric/MixedEquilibriumKineticReactions.hpp
     reactions/bulkGeneric/MixedEquilibriumKineticReactions_impl.hpp
     reactions/bulkGeneric/Parameters.hpp
     reactions/bulkGeneric/GenericSystems.hpp
     reactions/bulkGeneric/SpeciesUtilities.hpp
     reactions/unitTestUtilities/equilibriumReactionsTestUtilities.hpp
     reactions/unitTestUtilities/kineticReactionsTestUtilities.hpp
     reactions/unitTestUtilities/mixedReactionsTestUtilities.hpp
     reactions/massActions/MassActions.hpp
     )

set( hpcReact_sources)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)



#target_link_libraries(MyExecutable PRIVATE LAPACK::LAPACK)
     
set( hpcReack_dependencies
     LAPACK::LAPACK
     BLAS::BLAS
     )

blt_add_library( NAME             hpcReact
#                 SOURCES          ${hpcReact_sources}
                 HEADERS          ${hpcReact_headers}
                 DEPENDS_ON       ${hpcReact_dependencies}
                )

target_include_directories( hpcReact
                            INTERFACE
                            $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}../include>
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                            $<INSTALL_INTERFACE:include> )

blt_print_target_properties( TARGET hpcReact )
message(STATUS "HPCReact/src CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

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

add_subdirectory( reactions/bulkGeneric/unitTests )
add_subdirectory( reactions/geochemistry/unitTests )
add_subdirectory( reactions/massActions/unitTests )
add_subdirectory( common/unitTests )
add_subdirectory( docs )

hpcReact_add_code_checks( PREFIX hpcReact
                          UNCRUSTIFY_CFG_FILE ${PROJECT_SOURCE_DIR}/src/uncrustify.cfg
                          EXCLUDES "blt/*" )

