add_library(meshoui)

target_sources(meshoui PRIVATE
        mesh.cpp
        polygon.cpp
        plane.cpp
        clipping_surfaces.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
        )

target_include_directories(meshoui PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
target_link_libraries(meshoui PUBLIC OpenMeshCore MathUtils::MathUtils)
add_dependencies(meshoui check_git_${PROJECT_NAME}) # For git_watcher to fetch git informations before effective build

if (MESHOUI_USE_MMG)
    target_sources(meshoui PRIVATE remesher.cpp)
    if (BUILD_SHARED_LIBS)
        target_link_libraries(meshoui PRIVATE Mmg::libmmgs_so)
    else ()
        target_link_libraries(meshoui PRIVATE Mmg::libmmgs_a)
    endif ()
    target_compile_definitions(meshoui PUBLIC USE_MMG)
endif ()

if (USE_VTK)
    target_sources(meshoui PRIVATE vtkmesh.cpp)
    target_link_libraries(meshoui PRIVATE ${VTK_LIBRARIES})
    target_compile_definitions(meshoui PUBLIC USE_VTK)
endif ()
