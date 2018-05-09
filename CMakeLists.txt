# HelloWorld is a minimal sample creating, stepping and deleting a Bullet dynamics world

set (CMAKE_CXX_STANDARD 11)
cmake_minimum_required(VERSION 3.5)

INCLUDE_DIRECTORIES(
src
main
src/maths
)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# ADD_DEFINITION(-g -pg)

ADD_LIBRARY(
imrt
src/Station.cpp
src/EvaluationFunction.cpp
src/maths/Matrix.cpp
)


ADD_EXECUTABLE(TEST main/test.cpp)

target_link_libraries(TEST imrt)

IF (INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
			SET_TARGET_PROPERTIES(AppHelloWorld PROPERTIES  DEBUG_POSTFIX "_Debug")
			SET_TARGET_PROPERTIES(AppHelloWorld PROPERTIES  MINSIZEREL_POSTFIX "_MinsizeRel")
			SET_TARGET_PROPERTIES(AppHelloWorld PROPERTIES  RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
ENDIF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
