cmake_minimum_required(VERSION 3.13)
project(final_proj)

set(CMAKE_CXX_STANDARD 20)
add_definitions(-w)

find_package(LAPACK REQUIRED)
include_directories(/usr/include/)
link_libraries(/usr/lib/x86_64-linux-gnu/liblapacke.so)
link_libraries(${LAPACK_LIBRARIES})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})


find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("Found OpenMP")
    set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS})
endif()

add_executable(final_proj src/main.cpp src/Input.cpp src/Input.h)
target_link_libraries(final_proj  ${Boost_LIBRARIES})


