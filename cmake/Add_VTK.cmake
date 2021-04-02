message(STATUS "******* FINDING VTK dependency from ${PROJECT_NAME} (minimal requested version: ${VTK_TAG}) *******")
find_package(VTK ${VTK_TAG} REQUIRED NO_MODULE)
include(${VTK_USE_FILE})
