#
FILE(GLOB SOURCES *.f)

add_library(glimmersolveslapcmake ${SOURCES})
include_directories ( ./ )

