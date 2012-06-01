#
FILE(GLOB SOURCES *.F90)

add_library(glimmersolvecmake ${SOURCES})
add_dependencies(glimmersolvecmake glimmermodules)

include_directories (../ ${GLIMMER_BINARY_DIR}/libglimmer )

add_subdirectory(SLAP)
