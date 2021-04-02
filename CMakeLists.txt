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
set(CMAKE_CXX_STANDARD 17)

#-----------------------------------------------------------------------------
# Setting the configuration of the build (Debug, Release...)
#-----------------------------------------------------------------------------

# Access to local cmake scripts
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#=============================================================================
# MESHOUI options
#=============================================================================
option(MESHOUI_BUILD_TESTS "Activate build tests" ON)
option(MESHOUI_BUILD_TOOLS "Activate build tools" ON)
option(MESHOUI_USE_VTK "Enable VTK Support" ON)


#=============================================================================
# Retrieving the current Git revision
#=============================================================================
set(PRE_CONFIGURE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/meshoui/version.cpp.in CACHE FILEPATH "" FORCE)
set(POST_CONFIGURE_FILE ${CMAKE_CURRENT_BINARY_DIR}/src/meshoui/version.cpp CACHE FILEPATH "" FORCE)
set(GIT_WORKING_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE FILEPATH "" FORCE)
set(GIT_STATE_FILE ${CMAKE_CURRENT_BINARY_DIR}/git-state-hash CACHE FILEPATH "" FORCE)
include(cmake/git_watcher.cmake)

if (NOT EXISTS ${POST_CONFIGURE_FILE})
    configure_file(${PRE_CONFIGURE_FILE} ${POST_CONFIGURE_FILE})
endif()



# =============================================================================
# MESHOUI DEPENDENCIES
# =============================================================================

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/URL.conf.cmake)

include(Add_OpenMesh)
include(Add_MathUtils)
include(Add_mmg)


message(STATUS "MESHOUI USE VTK : ${MESHOUI_USE_VTK}" )

if (MESHOUI_USE_VTK)
    include(Add_VTK)
endif()

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

#=============================================================================
# Adding tools
#=============================================================================
if (${MESHOUI_BUILD_TOOLS})
    include(Add_cxxopt)
    add_subdirectory(tools)
endif()
