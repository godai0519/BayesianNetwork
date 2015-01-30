cmake_minimum_required(VERSION 2.8)

add_library(bp STATIC
  ../src/bp.cpp
)

add_library(graph STATIC
  ../src/graph.cpp
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

add_library(likelihood_weighting STATIC
  ../src/likelihood_weighting.cpp
)

add_library(utility STATIC
  ../src/utility.cpp
)

install(TARGETS bp DESTINATION lib)
install(TARGETS graph DESTINATION lib)
install(TARGETS cpt DESTINATION lib)
install(TARGETS matrix DESTINATION lib)
install(TARGETS sampling DESTINATION lib)
install(TARGETS likelihood_weighting DESTINATION lib)
install(TARGETS utility DESTINATION lib)
