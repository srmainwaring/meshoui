
if (MESHOUI_USE_MMG)
    add_executable(meshouix meshouix.cpp meshouix.cpp)
    target_link_libraries(meshouix meshoui cxxopts)
endif ()
