
# TODO : ajouter les fichiers source cpp dans la variable suivante :
set(MESHOUI_SRC_FILES
        mesh.cpp
        vtkmesh.cpp
        remesher.cpp
        )


add_library(meshoui SHARED)

target_sources(meshoui PRIVATE ${MESHOUI_SRC_FILES})
target_link_libraries(meshoui
        OpenMeshCore
        MathUtils::MathUtils
        ${VTK_LIBRARIES}
        Mmg::libmmgs_so
        )

target_include_directories(meshoui PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
