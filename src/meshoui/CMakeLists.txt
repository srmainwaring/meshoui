
# TODO : ajouter les fichiers source cpp dans la variable suivante :
set(MESHOUI_SRC_FILES
    loader.cpp
        )


add_library(meshoui STATIC)

target_sources(meshoui PRIVATE ${MESHOUI_SRC_FILES})
target_link_libraries(meshoui PRIVATE 
	)

target_include_directories(meshoui PUBLIC ${CMAKE_CURRENT_LIST_DIR}/..)
