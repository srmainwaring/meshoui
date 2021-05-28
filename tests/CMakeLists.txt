
message(STATUS "    ...meshoui TESTS :")

set(MESHOUI_TESTS
        test_analytical_surface_integration
        test_clipper
        test_load_mesh
        test_numerical_integration
        test_plane
        test_polygon
        test_symmetry_horizontal_plane
        test_merger
        )

foreach (test ${MESHOUI_TESTS})
    message(STATUS "        ...add ${test}")

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} PUBLIC meshoui gtest gtest_main)
    if (USE_VTK)
        target_link_libraries(${test} PRIVATE ${VTK_LIBRARIES})
    endif ()
    if (MESHOUI_USE_MMG)
        if (BUILD_SHARED_LIBS)
            target_link_libraries(${test} PRIVATE  Mmg::libmmgs_so)
        else ()
            target_link_libraries(${test} PRIVATE  Mmg::libmmgs_a)
        endif ()
    endif()

    add_test(NAME ${test} COMMAND ${test})

endforeach ()

# Data
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
