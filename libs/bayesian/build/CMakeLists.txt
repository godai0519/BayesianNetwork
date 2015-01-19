cmake_minimum_required(VERSION 2.8)

add_library(bp STATIC
  ../src/bp.cpp
)

add_library(cpt STATIC
  ../src/cpt.cpp
)

add_library(matrix STATIC
  ../src/matrix.cpp
)

add_library(sampling STATIC
  ../src/sampling.cpp
)

install(TARGETS bp DESTINATION lib)
install(TARGETS cpt DESTINATION lib)
install(TARGETS matrix DESTINATION lib)
install(TARGETS sampling DESTINATION lib)
