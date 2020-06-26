
message(STATUS "    ...meshoui TESTS :")

set(MESHOUI_TESTS
        test_load_mesh
        test_integration
        test_symmetry_horizontal_plane
        )

 foreach (test ${MESHOUI_TESTS})
    message(STATUS "        ...add ${test}")

    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} meshoui)
    if (MESHOUI_USE_VTK)
        target_compile_definitions(${test} PUBLIC -DMESHOUI_USE_VTK)
    endif ()

endforeach ()

set(subdirs
        # add here the test source subdirectories
        )

foreach (subdir ${subdirs})
    add_subdirectory(${subdir})
endforeach ()

## TODO : pour le moment pour chaque target, ajouter les 2 lignes suivantes
#
#add_executable(test_load_mesh test_load_mesh.cpp)
#target_link_libraries(test_load_mesh meshoui)
#
#add_executable(test_integration test_integration.cpp)
#target_link_libraries(test_integration meshoui)
#
#add_executable(test_symmetry_horizontal_plane test_symmetry_horizontal_plane.cpp)
#target_link_libraries(test_symmetry_horizontal_plane meshoui)