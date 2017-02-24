cmake_minimum_required(VERSION 2.8)
project(SmoothLifeTeachingExample CXX)

set (CMAKE_CXX_STANDARD 11)

include_directories(src)
include_directories(SYSTEM lib)

add_subdirectory(src)
add_subdirectory(test)

