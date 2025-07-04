#
set( PREPROCESSOR_DEFINES CUDA
                          HIP
                        )

set( USE_CONFIGFILE ON CACHE BOOL "" )
foreach( DEP in ${PREPROCESSOR_DEFINES})
    if( ${DEP}_FOUND OR ENABLE_${DEP} )
        set( SHIVA_USE_${DEP} TRUE )
    endif()
endforeach()



# configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/src/ShivaConfig.hpp.in
#                 ${CMAKE_BINARY_DIR}/include/ShivaConfig.hpp )

# Install the generated header.
# install( FILES ${CMAKE_BINARY_DIR}/include/ShivaConfig.hpp
#          DESTINATION include )

# Configure and install the CMake config
# configure_file( ${CMAKE_CURRENT_LIST_DIR}/shiva-config.cmake.in
#                 ${PROJECT_BINARY_DIR}/share/shiva/cmake/shiva-config.cmake)

# Set up cmake package config file

# set(SHIVA_INSTALL_INCLUDE_DIR "include" CACHE STRING "")
# set(SHIVA_INSTALL_CONFIG_DIR "lib" CACHE STRING "")
# set(SHIVA_INSTALL_LIB_DIR "lib" CACHE STRING "")
# set(SHIVA_INSTALL_BIN_DIR "bin" CACHE STRING "")
# set(SHIVA_INSTALL_CMAKE_MODULE_DIR "${SHIVA_INSTALL_CONFIG_DIR}/cmake" CACHE STRING "")
# set(SHIVA_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "" FORCE)


# include(CMakePackageConfigHelpers)
# configure_package_config_file(
#     ${CMAKE_CURRENT_SOURCE_DIR}/cmake/shiva-config.cmake.in
#     ${CMAKE_CURRENT_BINARY_DIR}/shiva-config.cmake
#   INSTALL_DESTINATION
#     ${SHIVA_INSTALL_CONFIG_DIR}
#   PATH_VARS
#     SHIVA_INSTALL_INCLUDE_DIR
#     SHIVA_INSTALL_LIB_DIR
#     SHIVA_INSTALL_BIN_DIR
#     SHIVA_INSTALL_CMAKE_MODULE_DIR
#   )


# install( FILES ${CMAKE_CURRENT_BINARY_DIR}/shiva-config.cmake
#          DESTINATION share/shiva/cmake/)
