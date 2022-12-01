cmake_minimum_required(VERSION 3.18)

project(ForceField)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(RDKit REQUIRED)
find_package(Boost COMPONENTS timer system REQUIRED)
add_executable(input input.cpp)
target_link_libraries(input RDKit::SmilesParse RDKit::FileParsers RDKit::GraphMol RDKit::RDGeneral
   Boost::timer)