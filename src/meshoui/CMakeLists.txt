
set(MESHOUI_SRC_FILES
        mesh.cpp
#        remesher.cpp
        polygon.cpp
        plane.cpp
        clipping_surfaces.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
        )


if (USE_VTK)
    list(APPEND MESHOUI_SRC_FILES vtkmesh.cpp)
endif ()

add_library(meshoui SHARED)
add_dependencies(meshoui check_git_${PROJECT_NAME}) # For git_watcher to fetch git informations before effective build

target_sources(meshoui PRIVATE ${MESHOUI_SRC_FILES})


target_link_libraries(meshoui PRIVATE
        MathUtils::MathUtils
        )
target_link_libraries(meshoui PUBLIC OpenMeshCore)

if (MESHOUI_USE_MMG)
    target_sources(meshoui PRIVATE remesher.cpp)
    target_link_libraries(meshoui PRIVATE Mmg::libmmgs_so)
    target_compile_definitions(meshoui PRIVATE USE_MMG)
endif ()

if (USE_VTK)
    target_link_libraries(meshoui PRIVATE ${VTK_LIBRARIES})
    target_compile_definitions(meshoui PRIVATE -DUSE_VTK)
endif ()

target_include_directories(meshoui PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
