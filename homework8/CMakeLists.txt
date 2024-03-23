cmake_minimum_required(VERSION 3.12)
project(HW8 VERSION 1.0 LANGUAGES C)

find_package(OpenMP)

add_executable(HW8 HW8.c)

include_directories(${MPI_INCLUDE_PATH}) # Addition specific to MacOS

target_link_libraries(HW8 PUBLIC
	OpenMP::OpenMP_C)
