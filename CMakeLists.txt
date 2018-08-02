# HelloWorld is a minimal sample creating, stepping and deleting a Bullet dynamics world

set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5)
project(imrt)


set(SOURCES src/solution/EvaluationFunction.cpp src/solution/Station.cpp src/base/Collimator.cpp src/tools/Matrix.cpp src/base/Volume.cpp src/solution/Plan.cpp src/ils/ApertureILS.cpp src/ils/IntensityILS.cpp)
set(HEADERS src/solution/EvaluationFunction.h src/solution/Plan.h src/solution/Station.h src/base/Collimator.h src/base/Volume.h src/ils/ILS.h src/ils/ApertureILS.h src/ils/IntensityILS.h)

INCLUDE_DIRECTORIES(
src
main
src/tools
src/base
src/solution
src/ils
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# ADD_DEFINITION(-g -pg)

ADD_LIBRARY(imrt ${SOURCES})


ADD_EXECUTABLE(TEST main/test.cpp)
ADD_EXECUTABLE(IAS main/intensity_aperture_solver.cpp)
#ADD_EXECUTABLE(AS main/aperture_solver.cpp)

target_link_libraries(TEST imrt)
target_link_libraries(IAS imrt)
#target_link_libraries(AS imrt)
