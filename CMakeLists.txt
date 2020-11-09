cmake_minimum_required(VERSION 3.16)

project(Root_simulation CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT)

include(${ROOT_USE_FILE})

add_compile_options(-Wall -Wextra)

add_subdirectory(generation)
add_subdirectory(analysis)
add_subdirectory(tex)