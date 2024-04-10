cmake_minimum_required(VERSION 3.12)
project(HW10 VERSION 1.0 LANGUAGES C)

find_package(OpenMP)

add_executable(HW10 HW10.c)

include_directories(${MPI_INCLUDE_PATH}) # Addition specific to MacOS

target_link_libraries(HW10 PUBLIC
	OpenMP::OpenMP_C
	m)
