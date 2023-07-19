cmake_minimum_required(VERSION 3.10)

project(cnvCompare VERSION 1.4.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -O3 -Wall -DBOOST_LOG_DYN_LINK "  CACHE STRING "Set C++ Compiler Flags" FORCE)

find_package(Boost 1.46.1 REQUIRED COMPONENTS program_options)
# check if boost was found
if(Boost_FOUND)
    message ("boost found")
else()
    message (FATAL_ERROR "Cannot find Boost")
endif()

find_package(Boost 1.46.1 REQUIRED COMPONENTS log)
# check if boost was found
if(Boost_FOUND)
    message ("boost found")
else()
    message (FATAL_ERROR "Cannot find Boost")
endif()

find_package(Boost 1.46.1 REQUIRED COMPONENTS log_setup)
# check if boost was found
if(Boost_FOUND)
    message ("boost found")
else()
    message (FATAL_ERROR "Cannot find Boost")
endif()

find_package(Boost 1.46.1 REQUIRED COMPONENTS filesystem)
# check if boost was found
if(Boost_FOUND)
    message ("boost found")
else()
    message (FATAL_ERROR "Cannot find Boost")
endif()

# target_include_directories("src")
add_executable(cnvCompare src/cnvCompare.cpp src/main.cpp src/utils.cpp)

# link against the boost libraries
target_link_libraries(cnvCompare PRIVATE Boost::program_options)
target_link_libraries(cnvCompare PRIVATE Boost::log)
target_link_libraries(cnvCompare PRIVATE Boost::log_setup)
target_link_libraries(cnvCompare PRIVATE Boost::filesystem)