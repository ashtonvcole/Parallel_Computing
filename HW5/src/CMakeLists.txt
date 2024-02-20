cmake_minimum_required(VERSION 3.12)
project(HW5 VERSION 1.0 LANGUAGES C)

find_package(MPI)

add_executable(HW5 HW5.c)

include_directories(${MPI_INCLUDE_PATH}) # Addition specific to MacOS

add_subdirectory(tools)

target_link_libraries(HW5 PUBLIC
	tools
	${MPI_C_LIBRARIES})
