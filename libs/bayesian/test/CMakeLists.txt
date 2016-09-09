add_executable(network_test network.cpp)

set_target_properties( network_test PROPERTIES OUTPUT_NAME "network" )

add_test( network_test network )
