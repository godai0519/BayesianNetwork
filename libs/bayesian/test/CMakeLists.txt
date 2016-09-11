include_directories( ${Boost_INCLUDE_DIRS} )
find_package( Boost COMPONENTS unit_test_framework REQUIRED )

add_executable(network_test network.cpp)
set_target_properties( network_test PROPERTIES OUTPUT_NAME "network" )
add_test( network_test network )
target_link_libraries(network_test ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

