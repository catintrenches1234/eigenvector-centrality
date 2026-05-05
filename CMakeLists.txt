cmake_minimum_required(VERSION 3.10)

project(eigenvector_centrality LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

add_executable(main
    src/main.cpp
    src/graph.cpp
    src/eigenvector_centrality.cpp
)

target_include_directories(main PRIVATE include)

target_link_libraries(main PRIVATE Eigen3::Eigen)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
        main PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -O2
    )
elseif (MSVC)
    target_compile_options(
        main PRIVATE
        /W4
        /O2
    )
endif()
