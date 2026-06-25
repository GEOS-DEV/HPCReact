# ============================================================================
# Auto-generate reaction system headers from JSON definitions
# ============================================================================
set(REACTION_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/reactions/data")
set(HPCREACT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(REACTION_GEN_SCRIPT "${HPCREACT_ROOT}/scripts/generate_reaction_system.py")
set(REACTION_VALIDATE_SCRIPT "${HPCREACT_ROOT}/scripts/validate_reaction_json.py")
set(GEOCHEMISTRY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/reactions/geochemistry")

# Find all JSON reaction definition files
file(GLOB REACTION_JSON_FILES "${REACTION_DATA_DIR}/*.json")

# Filter out example files (optional - remove this if you want to generate examples too)
# Commented out to generate all JSON files including examples
# list(FILTER REACTION_JSON_FILES EXCLUDE REGEX "example_.*\\.json$")

set(GENERATED_REACTION_HEADERS)

# Generate C++ header for each JSON file during CMake configuration
foreach(JSON_FILE ${REACTION_JSON_FILES})
    get_filename_component(BASENAME ${JSON_FILE} NAME_WE)

    # Convert filename to CamelCase for class name
    # e.g., my_system.json -> Mysystem.hpp
    string(SUBSTRING ${BASENAME} 0 1 FIRST_CHAR)
    string(TOUPPER ${FIRST_CHAR} FIRST_CHAR_UPPER)
    string(REGEX REPLACE "^.(.*)" "${FIRST_CHAR_UPPER}\\1" CLASSNAME ${BASENAME})

    set(HPP_FILE "${GEOCHEMISTRY_DIR}/${CLASSNAME}.hpp")

    # Step 1: Validate JSON first
    execute_process(
        COMMAND python3 ${REACTION_VALIDATE_SCRIPT} ${JSON_FILE}
        RESULT_VARIABLE VALIDATE_RESULT
        OUTPUT_QUIET
        ERROR_QUIET
    )

    if(NOT VALIDATE_RESULT EQUAL 0)
        message(WARNING "${BASENAME}.json validation failed - skipping generation")
        continue()
    endif()

    # Step 2: Generate C++ header (only if validation passed)
    message(STATUS "Generating reaction system: ${CLASSNAME}.hpp from ${BASENAME}.json")
    execute_process(
        COMMAND python3 ${REACTION_GEN_SCRIPT} ${JSON_FILE} -o ${HPP_FILE} --add-to-registry
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE GEN_RESULT
        OUTPUT_VARIABLE GEN_OUTPUT
        ERROR_VARIABLE GEN_ERROR
    )

    if(NOT GEN_RESULT EQUAL 0)
        message(WARNING "Failed to generate ${CLASSNAME}.hpp: ${GEN_ERROR}")
    else()
        message(STATUS "  -> Generated: ${HPP_FILE}")
        list(APPEND GENERATED_REACTION_HEADERS ${HPP_FILE})
    endif()
endforeach()

# Note: GeochemicalSystems.hpp is now updated incrementally via --add-to-registry flag above
# If you need to regenerate it from scratch, manually run:
#   python3 scripts/generate_reaction_system.py --generate-registry src/reactions
set(REGISTRY_HEADER "${GEOCHEMISTRY_DIR}/GeochemicalSystems.hpp")
if(EXISTS ${REGISTRY_HEADER})
    list(APPEND GENERATED_REACTION_HEADERS ${REGISTRY_HEADER})
endif()

# No need for custom target since files are generated during configuration
# The generated headers are already included in hpcReact_headers below

# ============================================================================

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
     ${GENERATED_REACTION_HEADERS}
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

if( NOT is_submodule )
  hpcReact_add_code_checks( PREFIX hpcReact
                            UNCRUSTIFY_CFG_FILE ${PROJECT_SOURCE_DIR}/src/uncrustify.cfg
                            EXCLUDES "blt/*" )
endif()
