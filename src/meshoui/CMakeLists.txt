
# TODO : ajouter les fichiers source cpp dans la variable suivante :
if (HDB5IO_USE_VTK)
	set(MESHOUI_SRC_FILES
		mesh.cpp
		vtkmesh.cpp
			)
	else()
	set(MESHOUI_SRC_FILES
			mesh.cpp
			)
endif()

add_library(meshoui SHARED)

target_sources(meshoui PRIVATE ${MESHOUI_SRC_FILES})
target_link_libraries(meshoui
	OpenMeshCore
	MathUtils::MathUtils
	)

if (MESHOUI_USE_VTK)
	target_link_libraries(meshoui ${VTK_LIBRARIES})
	target_compile_definitions(meshoui PUBLIC -DMESHOUI_USE_VTK)
endif ()

target_include_directories(meshoui PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
