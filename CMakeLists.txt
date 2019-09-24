#=============================================================================
#
# CMake configuration file for the MESHOUI project
#
#=============================================================================


#-----------------------------------------------------------------------------
# Check backward compatibility
#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.11)

#-----------------------------------------------------------------------------
# Prevent in-source builds
#-----------------------------------------------------------------------------
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message( FATAL_ERROR  "In-source build is not possible and not recommended. Choose an empty directory for build output.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})

#-----------------------------------------------------------------------------
# Define project name and specify version information
#-----------------------------------------------------------------------------
project(meshoui)
set(CMAKE_CXX_STANDARD 14)

#-----------------------------------------------------------------------------
# Setting the configuration of the build (Debug, Release...)
#-----------------------------------------------------------------------------

# Access to local cmake scripts
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#=============================================================================
# MESHOUI options
#=============================================================================
option(MESHOUI_BUILD_TESTS "Activate build tests" ON)


# =============================================================================
# MESHOUI DEPENDENCIES
# =============================================================================

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/URL.conf.cmake)

include(Add_OpenMesh)

#=============================================================================
# Adding the MESHOUI project sources
#=============================================================================
add_subdirectory(src/meshoui)

#=============================================================================
# Adding tests
#=============================================================================
if (${MESHOUI_BUILD_TESTS})
    enable_testing()
    include(Add_googletest)
    add_subdirectory(tests)
endif()
