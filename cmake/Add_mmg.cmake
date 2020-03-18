
include(FetchContent)

#    set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(mmg
        GIT_REPOSITORY ${mmg_URL}
        GIT_TAG ${mmg_TAG}
        )

FetchContent_GetProperties(mmg)
if(NOT mmg)
    message(STATUS "Downloading, Configuring and Generating 'mmg' dependency")
    FetchContent_Populate(mmg)

    # mmg BUILD OPTIONS
    set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

    add_subdirectory(${mmg_SOURCE_DIR} ${mmg_BINARY_DIR})

else()
    message(STATUS "mmg already populated")
endif()
