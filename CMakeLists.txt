# HelloWorld is a minimal sample creating, stepping and deleting a Bullet dynamics world

set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5)
project(imrt)


set(SOURCES src/EvaluationFunction.cpp src/Station.cpp src/maths/Matrix.cpp)
set(HEADERS src/EvaluationFunction.h src/Plan.h src/Station.h)

INCLUDE_DIRECTORIES(
src
main
src/maths
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# ADD_DEFINITION(-g -pg)

ADD_LIBRARY(imrt ${SOURCES})


ADD_EXECUTABLE(TEST main/test.cpp)

target_link_libraries(TEST imrt)

