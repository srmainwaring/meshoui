
find_package(OpenMesh QUIET)

if (NOT OpenMesh_FOUND)
    include(FetchContent)

    FetchContent_Declare(OpenMesh
            GIT_REPOSITORY ${OpenMesh_URL}
            GIT_TAG ${OpenMesh_TAG}
            )

    FetchContent_GetProperties(OpenMesh)
    if(NOT openmesh_POPULATED)
        message(STATUS "Downloading, Configuring and Generating 'OpenMesh' dependency")
        FetchContent_Populate(OpenMesh)

        # OpenMesh BUILD OPTIONS
        set(BUILD_APPS OFF CACHE BOOL "" FORCE)

        add_subdirectory(${openmesh_SOURCE_DIR} ${openmesh_BINARY_DIR})
    else()
        message(STATUS "OpenMesh already populated")
    endif()
endif()

if (TARGET OpenMeshCore)
    get_target_property(INC OpenMeshCore INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "Found OpenMesh : ${OPENMESH_INCLUDE_DIR}")
    target_include_directories(OpenMeshCore INTERFACE ${OPENMESH_INCLUDE_DIR})

else()
    message(STATUS "OpenMesh target NOT FOUND")
endif()
