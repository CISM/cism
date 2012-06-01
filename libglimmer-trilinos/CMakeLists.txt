#
FILE(GLOB SOURCES *.cpp)
FILE(GLOB HEADERS *.hpp)

add_library(glimmertrilinoscmake ${SOURCES} ${HEADERS})
include_directories ( ./ ../ ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

