################################################################
## Directories definition.
################################################################
set(SHARE_DIR "${CMAKE_SOURCE_DIR}/share")    # Directory containing test and configuration files.

##--------------- Do not edit below this line ----------------##

cmake_minimum_required(VERSION 2.8)
project(PROJECT_NAME CXX)

################################################################
## Set up some variables.
################################################################
set(MAIN "main")    # Name of the test source file to compile.

set(INCDIR  "${CMAKE_SOURCE_DIR}/include" )    # Directory where to find the library sources.
set(SRCDIR  "${CMAKE_SOURCE_DIR}/src" )    # Directory where to find the library sources.
set(TESTDIR "${CMAKE_SOURCE_DIR}/test")    # Directory where to find the test sources.
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

file(GLOB MAIN_SRC ${TESTDIR}/${MAIN}.cpp)
file(GLOB SRCS ${SRCDIR}/*.cpp)
file(GLOB HDRS ${INCDIR}/*.h)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -g")

################################################################
## Dependencies and third-party libraries.
################################################################
include_directories(SYSTEM "${INCDIR}")

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

################################################################
## Copy test and configuration files; create documentation dir.
################################################################
file(COPY ${SHARE_DIR} DESTINATION ${CMAKE_BINARY_DIR})

################################################################
## Target.
################################################################
add_executable(${MAIN} ${SRCS} ${MAIN_SRC})
