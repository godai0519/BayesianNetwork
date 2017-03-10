set(Boost_USE_STATIC_LIBS ON)
find_package( Boost COMPONENTS unit_test_framework REQUIRED )
include_directories( ${Boost_INCLUDE_DIRS} )

add_executable(network_test network.cpp)
set_target_properties( network_test PROPERTIES OUTPUT_NAME "network" )
add_test( network_test network )
target_link_libraries(network_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_executable(matrix_test matrix.cpp)
set_target_properties( matrix_test PROPERTIES OUTPUT_NAME "matrix" )
add_test( matrix_test matrix )
target_link_libraries(matrix_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
