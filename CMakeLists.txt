# HelloWorld is a minimal sample creating, stepping and deleting a Bullet dynamics world

set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5)
project(imrt)


set(SOURCES src/plan/EvaluationFunction.cpp src/plan/Station.cpp src/base/Collimator.cpp src/tools/Matrix.cpp src/base/Volume.cpp)
set(HEADERS src/plan/EvaluationFunction.h src/plan/Plan.h src/plan/Station.h src/base/Collimator.h src/base/Volume.h)

INCLUDE_DIRECTORIES(
src
main
src/tools
src/base
src/plan
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# ADD_DEFINITION(-g -pg)

ADD_LIBRARY(imrt ${SOURCES})


ADD_EXECUTABLE(TEST main/test.cpp)

target_link_libraries(TEST imrt)
