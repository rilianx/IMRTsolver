# HelloWorld is a minimal sample creating, stepping and deleting a Bullet dynamics world

set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5)
project(imrt)


set(SOURCES src/solution/EvaluationFunction.cpp src/solution/Station.cpp src/base/Collimator.cpp src/tools/Matrix.cpp src/base/Volume.cpp)
set(HEADERS src/solution/EvaluationFunction.h src/solution/Plan.h src/solution/Station.h src/base/Collimator.h src/base/Volume.h)

INCLUDE_DIRECTORIES(
src
main
src/tools
src/base
src/solution
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# ADD_DEFINITION(-g -pg)

ADD_LIBRARY(imrt ${SOURCES})


ADD_EXECUTABLE(TEST main/test.cpp)
ADD_EXECUTABLE(TEST2 main/test2.cpp)

target_link_libraries(TEST imrt)
target_link_libraries(TEST2 imrt)
