cmake_minimum_required(VERSION 2.8)

set(CMAKE_INSTALL_PREFIX "./")
set(CMAKE_CXX_FLAGS "-std=c++1y")

include_directories(./)
add_subdirectory(libs/bayesian/build)

