cmake_minimum_required(VERSION 2.8)

add_library(bp STATIC
  ../src/bp.cpp
)

add_library(matrix STATIC
  ../src/matrix.cpp
)

install(TARGETS bp DESTINATION lib)
install(TARGETS matrix DESTINATION lib)
