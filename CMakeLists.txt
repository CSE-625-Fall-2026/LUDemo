cmake_minimum_required(VERSION 3.16)
project(LUDemo VERSION 0.1.0 LANGUAGES CXX)

include(FetchContent)

FetchContent_Declare(matrix_class_demo
    GIT_REPOSITORY https://github.com/CSE-625-Fall-2026/MatrixClassDemo.git
    GIT_TAG v3.0.0
    GIT_SHALLOW TRUE
)
FetchContent_Declare(arbitrary_integer
    GIT_REPOSITORY https://github.com/AbstractClassroom-CPP/ArbitraryInteger.git
    GIT_TAG v4.0.0
    GIT_SHALLOW TRUE
)
FetchContent_Declare(rational_number
    GIT_REPOSITORY git@github.com:AbstractClassroom-CPP/RationalNumber.git
    GIT_TAG v5.0.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(matrix_class_demo arbitrary_integer rational_number)

add_executable(double_rref programs/double_rref.cpp)
target_link_libraries(double_rref PRIVATE MatrixClassDemo::matrix)

add_executable(double_lu programs/double_lu.cpp)
target_link_libraries(double_lu PRIVATE MatrixClassDemo::matrix)

add_executable(rational_rref programs/rational_rref.cpp)
target_link_libraries(rational_rref PRIVATE
    MatrixClassDemo::matrix
    ArbitraryInteger::ArbitraryInteger
    RationalNumber::RationalNumber
)

include(CTest)
if(BUILD_TESTING)
    add_test(NAME double_rref COMMAND double_rref)
    add_test(NAME double_lu COMMAND double_lu)
    add_test(NAME rational_rref COMMAND rational_rref)
endif()
