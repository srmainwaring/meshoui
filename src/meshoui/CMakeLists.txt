
set(MESHOUI_SRC_FILES
        mesh.cpp
        remesher.cpp
        polygon.cpp
        plane.cpp
        clipping_surfaces.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/version.cpp
        )

if (MESHOUI_USE_VTK)
    list(APPEND MESHOUI_SRC_FILES vtkmesh.cpp)
endif ()

add_library(meshoui SHARED)
add_dependencies(meshoui check_git_${PROJECT_NAME}) # For git_watcher to fetch git informations before effective build

target_sources(meshoui PRIVATE ${MESHOUI_SRC_FILES})
target_link_libraries(meshoui
        OpenMeshCore
        MathUtils::MathUtils
        Mmg::libmmgs_so
        )

if (MESHOUI_USE_VTK)
    target_link_libraries(meshoui ${VTK_LIBRARIES})
    target_compile_definitions(meshoui PUBLIC -DMESHOUI_USE_VTK)
endif ()

target_include_directories(meshoui PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
