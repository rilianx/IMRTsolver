# HelloWorld is a minimal sample creating, stepping and deleting a Bullet dynamics world

set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5)
project(imrt)


set(SOURCES src/solution/EvaluationFunction.cpp src/solution/Plan.cpp src/solution/Station.cpp
            src/base/Collimator.cpp src/tools/Matrix.cpp src/base/Volume.cpp src/ils/ApertureILS.cpp
            src/ils/IntensityILS2.cpp src/ils/IntensityGenerator.cpp)
set(HEADERS src/solution/EvaluationFunction.h src/solution/Plan.h src/solution/Station.h
            src/base/Collimator.h src/base/Volume.h src/ils/ILS.h src/ils/ApertureILS.h
            src/ils/IntensityILS2.h src/ils/IntensityGenerator.h)

INCLUDE_DIRECTORIES(
src
main
src/tools
src/base
src/solution
src/ils
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

# ADD_DEFINITION(-g -pg)

ADD_LIBRARY(imrt ${SOURCES})


#ADD_EXECUTABLE(TEST src/main/test.cpp)
ADD_EXECUTABLE(AS src/main/solver.cpp)
ADD_EXECUTABLE(DAO_ILS src/main/daoils.cpp)

#target_link_libraries(TEST imrt)
target_link_libraries(AS imrt)
target_link_libraries(DAO_ILS imrt)
