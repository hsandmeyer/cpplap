cmake_minimum_required(VERSION 3.0.0)

project("cpplap")

set(CMAKE_CXX_FLAGS "--std=c++17 -O3 -Wall -Wextra -Werror") 

include_directories("../")
add_executable(cpplap_test cpplap_test.cpp)
add_executable(example example.cpp)
