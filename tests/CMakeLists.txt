
set(TESTS
        test_load_mesh
        test_numerical_integration
        test_symmetry_horizontal_plane
        test_clipper
        test_mmg
        test_plane
        test_polygon
        )

foreach (test ${TESTS})
    message(STATUS "...add ${test}")

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} meshoui)
endforeach ()

# Data
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
