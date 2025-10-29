
set( hpcReact_headers
     common/macros.hpp
     common/CArrayWrapper.hpp
     reactions/exampleSystems/BulkGeneric.hpp
     reactions/geochemistry/Carbonate.hpp
     reactions/geochemistry/Forge.hpp
     reactions/geochemistry/GeochemicalSystems.hpp
     reactions/geochemistry/Ultramafics.hpp
     reactions/massActions/MassActions.hpp
     reactions/reactionsSystems/EquilibriumReactions.hpp
     reactions/reactionsSystems/EquilibriumReactionsAggregatePrimaryConcentration_impl.hpp
     reactions/reactionsSystems/EquilibriumReactionsReactionExtents_impl.hpp
     reactions/reactionsSystems/KineticReactions.hpp
     reactions/reactionsSystems/KineticReactions_impl.hpp
     reactions/reactionsSystems/MixedEquilibriumKineticReactions.hpp
     reactions/reactionsSystems/MixedEquilibriumKineticReactions_impl.hpp
     reactions/reactionsSystems/Parameters.hpp
     reactions/unitTestUtilities/equilibriumReactionsTestUtilities.hpp
     reactions/unitTestUtilities/kineticReactionsTestUtilities.hpp
     reactions/unitTestUtilities/mixedReactionsTestUtilities.hpp
     )

set( hpcReact_sources)

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

     
set( hpcReact_dependencies
     LAPACK::LAPACK
     BLAS::BLAS
     )

if( ENABLE_CUDA )
  list( APPEND hpcReact_dependencies cuda )
endif()

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

add_subdirectory( reactions/exampleSystems/unitTests )
add_subdirectory( reactions/geochemistry/unitTests )
add_subdirectory( reactions/massActions/unitTests )
add_subdirectory( common/unitTests )
add_subdirectory( docs )

hpcReact_add_code_checks( PREFIX hpcReact
                          UNCRUSTIFY_CFG_FILE ${PROJECT_SOURCE_DIR}/src/uncrustify.cfg
                          EXCLUDES "blt/*" )

