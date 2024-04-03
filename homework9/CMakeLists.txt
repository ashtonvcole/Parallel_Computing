set(CMAKE_C_COMPILER gcc)

cmake_minimum_required(VERSION 3.12)
project(HW9 VERSION 1.0 LANGUAGES C)

find_package(OpenMP)

add_executable(HW7 HW7.c)
add_executable(HW9 HW9.c)

include_directories(${MPI_INCLUDE_PATH}) # Addition specific to MacOS

target_link_libraries(HW7 PUBLIC
	OpenMP::OpenMP_C
	m)
target_link_libraries(HW9 PUBLIC
	OpenMP::OpenMP_C
	m)
