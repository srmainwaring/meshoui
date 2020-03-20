

# TODO : pour le moment pour chaque target, ajouter les 2 lignes suivantes

set(TESTS
        test_load_mesh
        test_integration
        test_symmetry_horizontal_plane
        test_clipper
        test_mmg
        )

foreach (test ${TESTS})
    message(STATUS "...add ${test}")

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} meshoui)
endforeach ()

# Data
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


#add_executable(test_load_mesh test_load_mesh.cpp)
#target_link_libraries(test_load_mesh meshoui)
#
#add_executable(test_integration test_integration.cpp)
#target_link_libraries(test_integration meshoui)
#
#add_executable(test_symmetry_horizontal_plane test_symmetry_horizontal_plane.cpp)
#target_link_libraries(test_symmetry_horizontal_plane meshoui)
