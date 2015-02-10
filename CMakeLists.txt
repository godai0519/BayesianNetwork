cmake_minimum_required(VERSION 2.8)
project(bayesian_network)
enable_testing()

find_package( Boost REQUIRED )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${Boost_INCLUDE_DIRS} )

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -O2")
set(CMAKE_C_FLAGS_RELEASE "-Wall -Wextra -g")

if( MSVC OR MSVC_IDE )
    if( MSVC_VERSION LESS 1800 )
        message(FATAL_ERROR "The library requires C++11 or more."
                            "Need at least Visual Studio 12.0 (a.k.a. Visual Studio 2013) "
                            "with Microsoft Visual C++ Compiler November 2013 CTP (v120_CTP_Nov2013)" )
    elseif( MSVC_VERSION EQUAL 1800 )
        message( "VC12: Use Microsoft Visual Studio 2013 with Microsoft Visual C++ Compiler November 2013 CTP (v120_CTP_Nov2013)" )
        set(CMAKE_GENERATOR_TOOLSET "CTP_Nov2013" CACHE STRING "Platform Toolset" FORCE)
    else()
    endif()
    
    add_compile_options( -DNOMINMAX )
else()
    add_compile_options( -std=c++1y )
endif()
subdirs( bayesian libs )
